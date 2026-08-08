//------------------------------------------------------------------------------------
// Copyright (C) 2003  Columbia Genome Center * All Rights Reserved. *
//
// Modifications by S.V. Rice, 2017
//------------------------------------------------------------------------------------

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <functional>
#include <ios>
#include <iterator>
#include <new>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include "matrix.h"
#include "util.h"

int maxNpar = 0; // maximum observed value of npar

//------------------------------------------------------------------------------------

class ArrayValuePair
{
public:
   ArrayValuePair(int id, double v)
      : arrayId(id), value(v) { }

   int arrayId;  // getId(), set()
   double value; // getValue(), set()
};

typedef std::binary_function<ArrayValuePair, ArrayValuePair, bool>
   ArrayValuePairBinaryFunction;

class SortIncreasing_ArrayValuePair : ArrayValuePairBinaryFunction
{
public:
   bool operator()(const ArrayValuePair& a, const ArrayValuePair& b) const
   {
      return (a.value < b.value);
   }
};

class SortDecreasing_ArrayValuePair : ArrayValuePairBinaryFunction
{
public:
   bool operator()(const ArrayValuePair& a, const ArrayValuePair& b) const
   {
      return (a.value > b.value);
   }
};

//------------------------------------------------------------------------------------

class RankedObservation
{
public:
   RankedObservation()
      : value(), position() { }

   RankedObservation(double inValue, int inPosition)
      : value(inValue), position(inPosition) { }

   double value;
   int position;
};

class Sort_RankedObservation
{
public:
   bool operator()(const RankedObservation& a, const RankedObservation& b) const
   {
      return (a.value < b.value ||
              (a.value == b.value && a.position < b.position));
   }
};

//------------------------------------------------------------------------------------

// One workspace belongs to one createEdgeMatrix() invocation.  Every buffer range
// that can be read must be reset or overwritten before each MI pair or partition.
class AdaptivePartitionWorkspace
{
public:
   AdaptivePartitionWorkspace()
      : observationCount(0), partitionLimit(0), poc(), kon(), poradi(), marg(),
        apor(), quadrant(), NN(), amarg() { }

   void initialize(int N, int M)
   {
      observationCount = N;
      partitionLimit   = M;

      poc.resize(M);
      kon.resize(M);
      poradi.resize(N);
      marg.resize(4 * static_cast<std::size_t>(M));
      apor.resize(N);
      quadrant.resize(N);
      NN.resize(4);
      amarg.resize(16);
   }

   bool matches(int N, int M) const
   {
      return observationCount == N && partitionLimit == M;
   }

   void resetEdge()
   {
      std::fill(poc.begin(), poc.end(), 1);
      std::fill(kon.begin(), kon.end(), observationCount);
      std::fill(marg.begin(), marg.end(), 0);

      for (int i = 0; i < observationCount; i++)
         poradi[i] = i + 1;

      marg[0] = marg[partitionLimit] = 1;
      marg[2 * partitionLimit] = marg[3 * partitionLimit] = observationCount;
   }

   int observationCount;
   int partitionLimit;
   std::vector<int> poc;
   std::vector<int> kon;
   std::vector<int> poradi;
   std::vector<int> marg;
   std::vector<int> apor;
   std::vector<unsigned char> quadrant;
   std::vector<int> NN;
   std::vector<int> amarg;
};

//------------------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const Node& n)
{
   out << "(" << n.mutinfo << ", " << n.intermediate << ")";
   return out;
}

//------------------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const Marker& m)
{
   out << "(" << m.label << "; " << m.accnum << "; " << m.idnum << "; "
       << (m.isActive ? "Active; " : "Not Active; ")
       << (m.isControl ? "Control" : "Not Control" ) << ")";
   return out;
}

//------------------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const Probe& p)
{
   out << "(" << p.value << ", " << p.pvalue << ")";
   return out;
}

//------------------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const Microarray_Set& ms)
{
   const char *delim = "\t";

   int numMarkers     = ms.markerset.size();
   int numMicroarrays = ms.uarrays.size();

   if (ms.header.size() > 0)
      std::copy(ms.header.begin(), ms.header.end(),
                std::ostream_iterator<std::string>(out, delim));
   else
   {
      // make up a col header
      out << "Id" << delim << "Desc" << delim;

      int col = 0;
      while (col < numMarkers - 1)
         out << "exp" << col++ << delim;

      out << "exp" << col;
   }

   out << "\n";

   for (int i = 0; i < numMarkers; i++)
   {
      out << ms.markerset[i].accnum << "\t" << ms.markerset[i].label;

      for (int j = 0; j < numMicroarrays; j++)
         out << delim << ms.uarrays[j][i].value;

      if (i < numMarkers - 1)
         out << "\n";
   }

   return out;
}

//------------------------------------------------------------------------------------

void Matrix::saveNode(int i, int j, double mi)
{
   if (!std::isfinite(mi))
      throw std::string("Refusing to store a non-finite MI value.");

   NodeMap& nmap = nmv[i];

   NodeMap::iterator npos = nmap.find(j);

   if (npos == nmap.end())
      nmap.insert(std::make_pair(j, Node(mi)));
   else
      npos->second.mutinfo = mi;
}

//------------------------------------------------------------------------------------

void Matrix::read(Microarray_Set& data, const Parameter& p)
{
   std::ifstream in(p.adjfile.c_str());
   if (!in.is_open())
      throw "Unable to open " + p.adjfile;

   read(in, data, p);

   in.close();
}

//------------------------------------------------------------------------------------

static double parseMutualInformation(const std::string& text,
                                     const std::string& source,
                                     const std::string& target,
                                     const std::string& filename)
{
   const char *begin = text.c_str();
   char *end = NULL;
   errno = 0;

   double mi = std::strtod(begin, &end);
   const bool converted = (end != NULL && end != begin);

   while (end != NULL && *end != '\0' &&
          std::isspace(static_cast<unsigned char>(*end)))
      end++;

   std::ostringstream s;

   if (!converted || end == NULL || *end != '\0' || errno == ERANGE)
   {
      s << "Invalid MI value \"" << text << "\" for edge " << source
        << " -> " << target << " in adjacency file \"" << filename << "\".";
      throw s.str();
   }

   if (!std::isfinite(mi))
   {
      s << "Non-finite MI value \"" << text << "\" for edge " << source
        << " -> " << target << " in adjacency file \"" << filename << "\".";
      throw s.str();
   }

   return mi;
}

//------------------------------------------------------------------------------------

void Matrix::read(std::istream& in, Microarray_Set& data, const Parameter& p)
{
   // The expression matrix defines the complete gene-ID space.  Pre-sizing
   // prevents sparse -j inputs from leaving requested rows out of bounds.
   nmv.assign(data.markerset.size(), NodeMap());
   adjacencyRowsPresent.assign(data.markerset.size(), false);

   std::string line;

   std::getline(in, line);
   while (line.length() > 0 && line[0] == '>')
      std::getline(in, line);

   while (in.good())
   {
      std::istringstream sin(line);
      std::string label;

      std::getline(sin, label, '\t');
      label = "_" + label;

      int geneId1 = data.getProbeId(label);
      if (geneId1 == -1)
         throw "Cannot find marker: " + label + " in the ADJ file!";

      data.markerset[geneId1].isActive = true;
      adjacencyRowsPresent[geneId1] = true;

      std::getline(sin, label, '\t');
      label = "_" + label;

      while (sin.good())
      {
         std::string value;

         std::getline(sin, value, '\t');

         const std::string source = data.markerset[geneId1].accnum.substr(1);
         const std::string target = label.substr(1);
         double mi = parseMutualInformation(value, source, target, p.adjfile);

         if (mi >= p.threshold)
         {
            int geneId2 = data.getProbeId(label);
            if (geneId2 == -1)
               throw "Cannot find marker: " + label + " in the ADJ file!";

            saveNode(geneId1, geneId2, mi);
         }

         std::getline(sin, label, '\t');
         label = "_" + label;
      }

      std::getline(in, line);
   }
}

//------------------------------------------------------------------------------------

bool Matrix::hasAdjacencyRow(int geneId) const
{
   return geneId >= 0 &&
          static_cast<std::size_t>(geneId) < adjacencyRowsPresent.size() &&
          adjacencyRowsPresent[geneId];
}

//------------------------------------------------------------------------------------

void Matrix::writeGeneLine(std::ostream& out, const Microarray_Set& data, int geneId)
{
   NodeMap& nmap = nmv[geneId];

   const Marker& marker = data.markerset[geneId];

   if (nmap.size() == 0 &&
       (!writeEmptyGenes || !marker.isActive && !marker.isControl))
      return;

   out << marker.accnum.substr(1);

   for (NodeMap::iterator npos = nmap.begin(); npos != nmap.end(); ++npos)
   {
      int   id   = npos->first;
      Node& node = npos->second;

      if (writeTriangular && id <= geneId || !writeReduced && node.intermediate >= 0)
         continue;

      const Marker& marker = data.markerset[id];

      out << "\t" << marker.accnum.substr(1);

      if (node.intermediate >= 0)
         out << "." << node.intermediate;

      out << "\t" << node.mutinfo;
   }

   out << std::endl;
}

//------------------------------------------------------------------------------------

void Matrix::writeGeneList(const Microarray_Set& data, const std::string& name,
                           int probeId)
{
   std::string filename = name + ".adj";

   std::ofstream out(filename.c_str());
   if (!out.is_open())
      throw "Unable to open " + filename;

   std::cout << "Writing gene list: "<< filename << std::endl;

   NodeMap& nmap = nmv[probeId];

   const Marker& marker = data.markerset[probeId];

   if (nmap.size() == 0 &&
       (!writeEmptyGenes || !marker.isActive && !marker.isControl))
      return;

   for (NodeMap::iterator npos = nmap.begin(); npos != nmap.end(); ++npos)
   {
      int   id   = npos->first;
      Node& node = npos->second;

      const Marker& marker = data.markerset[id];

      if (id != probeId)
         out << id << "\t" << marker.accnum << "t" << node.mutinfo << std::endl;
   }

   out.close();
}

//------------------------------------------------------------------------------------

void Matrix::write(std::ostream& out, const Microarray_Set& data,
                   const std::vector<int>& ids)
{
   int numIds = ids.size();

   if (numIds == 0)
   {
      int numNodeMaps = nmv.size();

      for (int id = 0; id < numNodeMaps; id++)
         writeGeneLine(out, data, id);
   }
   else
      for (int i = 0; i < numIds; i++)
         writeGeneLine(out, data, ids[i]);
}

//------------------------------------------------------------------------------------

void Matrix::write(const Microarray_Set& data, const std::vector<int>& ids,
                   const Parameter& p, bool writeFull)
{
   if (!writeFull)
      return;

   std::ofstream out(p.outfile.c_str());
   if (!out.is_open())
      throw "Unable to open " + p.outfile;

   std::cout << "Writing matrix: " << p.outfile << std::endl;

   out << ">  Input file      " << p.infile     << std::endl;
   //out << ">  ADJ file        " << p.adjfile    << std::endl;
   //out << ">  Output file     " << p.outfile    << std::endl;
   //out << ">  Kernel width    " << p.sigma      << std::endl;
   out << ">  MI threshold    " << p.threshold  << std::endl;
   out << ">  MI P-value      " << p.pvalue     << std::endl;
   out << ">  DPI tolerance   " << p.eps        << std::endl;
   //out << ">  Correction      " << p.correction << std::endl;
   out << ">  Subnetwork file " << p.subnetfile << std::endl;
   //out << ">  Hub probe       " << p.hub        << std::endl;
   //out << ">  Control probe   " << p.controlId  << std::endl;
   //out << ">  Condition       " << p.condition  << std::endl;
   //out << ">  Percentage      " << p.percent    << std::endl;
   //out << ">  TF annotation   " << p.annotfile  << std::endl;
   //out << ">  Filter mean     " << p.mean       << std::endl;
   //out << ">  Filter CV       " << p.cv         << std::endl;

   write(out, data, ids);

   out.close();

   std::cout << "Maximum observed npar: " << maxNpar << std::endl;
}

//------------------------------------------------------------------------------------

void Matrix::createEntries(int numEntries)
{
   // create entries for all rows in the matrix

   for (int i = 0; i < numEntries; i++)
      nmv.push_back(NodeMap());
}

//------------------------------------------------------------------------------------

void Matrix::addNode(int i, int j, double edgeValue, bool symmetric)
{
   saveNode(i, j, edgeValue);

   if (symmetric)
      saveNode(j, i, edgeValue);
}

//------------------------------------------------------------------------------------

double Matrix::getNodeMI(int geneId1, int geneId2)
{
   // this function returns a positive MI value if the node has been computed and
   // survived the thresholding; returns 0.0 if the node has been computed but did not
   // survive the thresholding; or returns -1.0 if the node has not yet been computed

   if (nmv.size() > geneId1 && nmv[geneId1].size() > 0)
   {
      NodeMap& nmap = nmv[geneId1];

      NodeMap::iterator npos = nmap.find(geneId2);

      return (npos == nmap.end() ? 0.0 : npos->second.mutinfo);
   }
   else if (nmv.size() > geneId2 && nmv[geneId2].size() > 0)
   {
      NodeMap& nmap = nmv[geneId2];

      NodeMap::iterator npos = nmap.find(geneId1);

      return (npos == nmap.end() ? 0.0 : npos->second.mutinfo);
   }
   else
      return -1.0;
}

//------------------------------------------------------------------------------------

static bool protectedByTFLogic(Transfac &transfac,
                               int geneId1, int geneId2, int geneId3)
{
   bool isA = (transfac.find(geneId1) != transfac.end());
   bool isB = (transfac.find(geneId2) != transfac.end());
   bool isC = (transfac.find(geneId3) != transfac.end());

   return ((isA || isB) && !isC);
}

//------------------------------------------------------------------------------------

void Matrix::reduceOneNode(int row_idx, double epsilon, Transfac& transfac)
{
   NodeMap& nmap = nmv[row_idx];

   std::vector<ArrayValuePair> miVector;

   for (NodeMap::iterator npos = nmap.begin(); npos != nmap.end(); ++npos)
      miVector.push_back(ArrayValuePair(npos->first, npos->second.mutinfo));

   SortDecreasing_ArrayValuePair sorter;
   std::sort(miVector.begin(), miVector.end(), sorter);

   int numPairs = miVector.size();

   for (int i = 0; i < numPairs; i++)
   {
      int    geneId1 = miVector[i].arrayId;
      double valueAB = miVector[i].value;

      double minMI   = valueAB / (1.0 - epsilon);

      for (int j = 0; j < i; j++)
      {
         int    geneId2 = miVector[j].arrayId;
         double valueAC = miVector[j].value;

         if (valueAC <= minMI)
            break;

         double valueBC = getNodeMI(geneId1, geneId2);

         if (valueBC > minMI && (transfac.size() == 0 ||
             !protectedByTFLogic(transfac, row_idx, geneId1, geneId2)))
         {
            NodeMap::iterator npos = nmap.find(geneId1);

            if (npos != nmap.end())
               npos->second.intermediate = geneId2;

            break;
         }
      }
   }
}

//------------------------------------------------------------------------------------

void Matrix::reduce(double epsilon, const std::vector<int>& ids, Transfac& transfac)
{
   std::time_t t1, t2;
   std::time(&t1);

   if (ids.size() == 0)
   {
      int numMarkers = nmv.size();

      for (int i = 0; i < numMarkers; i++)
         reduceOneNode(i, epsilon, transfac);
   }
   else
   {
      int numIds = ids.size();

      for (int i = 0; i < numIds; i++)
         reduceOneNode(ids[i], epsilon, transfac);
   }

   std::time(&t2);
   std::cout << "DPI running time is: " << std::difftime(t2, t1) << "\n";
}

//------------------------------------------------------------------------------------

int Microarray_Set::Get_Num_Active_Markers() const
{
   int numMarkers = markerset.size();
   int numActive  = 0;

   for (int i = 0; i < numMarkers; i++)
      if (markerset[i].isActive)
         numActive++;

   return numActive;
}

//------------------------------------------------------------------------------------

bool Microarray_Set::isSameGene(int i, int j) const
{
   return (i == j || markerset[i].label == markerset[j].label &&
                     markerset[i].label != "---");
}

//------------------------------------------------------------------------------------

int Microarray_Set::getAccessionId(const std::string& accnum) const
{
   if (!accessionIdsCurrent)
      rebuildAccessionIndex();

   std::unordered_map<std::string, int>::const_iterator found =
      accessionIds.find(accnum);

   if (found != accessionIds.end())
      return found->second;

   return -1; // did not find it
}

//------------------------------------------------------------------------------------

void Microarray_Set::rebuildAccessionIndex() const
{
   accessionIds.clear();
   accessionIds.reserve(markerset.size());

   int numMarkers = markerset.size();

   for (int i = 0; i < numMarkers; i++)
      // emplace() does not replace an existing entry. This preserves the legacy
      // linear lookup behavior when an expression matrix contains duplicate
      // accession IDs: the first expression row remains authoritative.
      accessionIds.emplace(markerset[i].accnum, i);

   accessionIdsCurrent = true;
}

//------------------------------------------------------------------------------------

int Microarray_Set::getProbeId(const std::string& label) const
{
   int len = label.length();

   for (int i = 0; i < len; i++)
      if (!std::isdigit(label[i]))     // label is not a number
         return getAccessionId(label); // perhaps it is an accession ID

   return std::atoi(label.c_str());    // label is a number, convert it to int
}

//------------------------------------------------------------------------------------

void Microarray_Set::Set_ColHeader(int i, const std::string& hdr)
{
   while (i >= header.size())
      header.push_back("");

   header[i] = hdr;
}

//------------------------------------------------------------------------------------

void Microarray_Set::Set_Marker(int i, const Marker& m)
{
   while (i >= markerset.size())
      markerset.push_back(Marker());

   markerset[i] = m;
   accessionIdsCurrent = false;
}

//------------------------------------------------------------------------------------

void Microarray_Set::Set_Probe(int i, int j, Probe p)
{
   while (i >= uarrays.size())
      uarrays.push_back(Microarray());

   while (j >= uarrays[i].size())
      uarrays[i].push_back(Probe());

   uarrays[i][j] = p;
}

//------------------------------------------------------------------------------------
// Read a matrix data file into a Microarray_Set. The file consists of sets of lines
// that start with a description of a marker (including accession number and
// descriptive label). The remainder of the lines are probe values (the values only,
// or (value, pvalue) pairs) for each microarray for the marker that started the line.
// We assign each marker an ID number according to the line number it appeared on.
// Conceptually, the file can be thought of as a matrix of probe data in column major
// order.
//
// The program can automatically tell whether the input file consists of
// (value, pvalue) pairs or single expression values. It will read the data and set
// up the structure correctly according to the format of the file. (If the data
// contains only expression measurements, a p-value of 0.0 will be assigned to each
// expression value.)
//------------------------------------------------------------------------------------

int Microarray_Set::readMarkerWithPvalue(std::istream& in, const int arrno)
{
   // Read a marker line with (value, pvalue) format

   int markern = 0; // number of expression values of each probe

   in.exceptions(std::ios_base::badbit | std::ios_base::failbit);

   try
   {
      do
      {
         std::string sPVal;
         double val, pval;

         in >> val >> sPVal;

         switch (sPVal[0])
         {
            case 'A': pval = 0.7; break;
            case 'M': pval = 0.5; break;
            case 'P': pval = 0.1; break;
            default : pval = std::atof(sPVal.c_str());
         }

         if (!std::isfinite(val))
         {
            std::ostringstream s;
            s << "Non-finite expression value at gene row " << arrno + 1
              << ", observation " << markern + 1 << ".";
            throw s.str();
         }

         if (!std::isfinite(pval))
         {
            std::ostringstream s;
            s << "Non-finite expression p-value at gene row " << arrno + 1
              << ", observation " << markern + 1 << ".";
            throw s.str();
         }

         Set_Probe(markern++, arrno, Probe(val, pval));
      }
      while (in.good() && in.get() != '\015' && in.peek() != EOF &&
             in.peek() != '\015');
   }
   catch (std::ios_base::failure& f)
   {
      std::ostringstream s;
      s << "Could not read data at line no: " << arrno + 2;
      std::cout << s.str() << std::endl;
      throw s.str();
   }

   return markern;
}

//------------------------------------------------------------------------------------

int Microarray_Set::readMarkerNoPvalue(std::istream& in, const int arrno)
{
   // Read a marker line with no pvalues

   int markern = 0; // number of expression values of each probe

   in.exceptions(std::ios_base::badbit | std::ios_base::failbit);

   try
   {
      do
      {
         double val;

         in >> val;

         if (!std::isfinite(val))
         {
            std::ostringstream s;
            s << "Non-finite expression value at gene row " << arrno + 1
              << ", observation " << markern + 1 << ".";
            throw s.str();
         }

         Set_Probe(markern++, arrno, Probe(val, 0.0));
      }
      while (in.good() && in.get() != '\015' && in.peek() != EOF &&
             in.peek() != '\015');
   }
   catch (std::ios_base::failure& f)
   {
      std::ostringstream s;
      s << "Could not read data at line no: " << arrno;
      std::cout << s.str() << std::endl;
      throw s.str();
   }

   return markern;
}

//------------------------------------------------------------------------------------

int Microarray_Set::readHeader(std::istream& in)
{
   in.exceptions(std::ios_base::badbit | std::ios_base::failbit);

   try
   {
      do
      {
         std::string hdr;

         std::getline(in, hdr, '\t');

         Set_ColHeader(header.size(), hdr);
      }
      while (in.good() && in.peek() != '\015' && in.peek() != EOF);
   }
   catch (std::ios_base::failure& f)
   {
      throw std::string("Error while reading file headers"
                        "(win/*nux/mac end of line?).");
   }

   return header.size() - 2;
}

//------------------------------------------------------------------------------------

static void validateObservationCount(int expected, int actual, int lineNumber)
{
   if (actual != expected)
   {
      std::ostringstream s;
      s << "Incorrect expression dimensions at line " << lineNumber
        << ": expected " << expected << " observations, found " << actual << ".";
      throw s.str();
   }
}

//------------------------------------------------------------------------------------

void Microarray_Set::read(const std::string& filename)
{
   std::ifstream in(filename.c_str());
   if (!in.is_open())
      throw "Unable to open " + filename;

   read(in);

   in.close();
}

//------------------------------------------------------------------------------------

void Microarray_Set::read(std::istream& in)
{
   std::ios_base::iostate oldState = in.exceptions();
   in.exceptions(std::ios_base::badbit | std::ios_base::failbit);

   try
   {
      std::string line;
      std::getline(in, line);
      std::istringstream sin(line);
      int arrno = readHeader(sin);

      if (arrno < 2)
      {
         std::ostringstream s;
         s << "Expression matrix must contain at least 2 observation columns; found "
           << arrno << ".";
         throw s.str();
      }

      // we need to decide whether the input file contain only expression
      // or (value, pvalue) pairs

      int bypass_line_cnt = 0;
      std::getline(in, line);

      while (line.length() >= 11 && line.substr(0, 11) == "Description")
      {
         bypass_line_cnt++;
         std::getline(in, line);
      }

      std::cout << "\n[READ] " << bypass_line_cnt
                << " Description lines bypassed." << std::endl;

      const int firstDataLine = bypass_line_cnt + 2;

      std::istringstream pin(line);
      std::vector<std::string> firstprobe;
      std::string token;

      do
      {
         std::getline(pin, token, '\t');
         firstprobe.push_back(token);
      }
      while (pin.good() && pin.peek() != '\015' && pin.peek() != EOF);

      int valueNo = firstprobe.size() - 2;

      if (arrno == valueNo)
      {
         // single expression values
         std::cout << "[READ] P-value columns not found." << std::endl;

         int proben = 0; // probe number

         std::istringstream fin(line);
         std::string accnum, label;

         std::getline(fin, accnum, '\t');
         std::getline(fin, label,  '\t');

         accnum = "_" + accnum;

         Marker m(proben, accnum, label);
         Set_Marker(proben, m);
         int observations = readMarkerNoPvalue(fin, proben);
         validateObservationCount(arrno, observations, firstDataLine);
         proben++;

         while (in.good() && in.peek() != EOF)
         {
            std::getline(in, line);
            std::istringstream sin(line);

            std::getline(sin, accnum, '\t');
            std::getline(sin, label,  '\t');

            accnum = "_" + accnum;

            Marker m(proben, accnum, label);
            Set_Marker(proben, m);

            int observations = readMarkerNoPvalue(sin, proben);
            validateObservationCount(arrno, observations,
                                     firstDataLine + proben);

            proben++;
         }
      }
      else if (2 * arrno == valueNo)
      {
         // (value, pvalue) pairs
         std::cout << "[READ] (value, p-value) pairs found." << std::endl;

         int proben = 0; // probe number

         std::istringstream fin(line);
         std::string accnum, label;

         std::getline(fin, accnum, '\t');
         std::getline(fin, label,  '\t');

         accnum = "_" + accnum;

         Marker m(proben, accnum, label);
         Set_Marker(proben, m);
         int observations = readMarkerWithPvalue(fin, proben);
         validateObservationCount(arrno, observations, firstDataLine);
         proben++;

         while (in.good() && in.peek() != EOF)
         {
            std::getline(in, line);
            std::istringstream sin(line);

            std::getline(sin, accnum, '\t');
            std::getline(sin, label,  '\t');

            accnum = "_" + accnum;

            Marker m(proben, accnum, label);
            Set_Marker(proben, m);

            int observations = readMarkerWithPvalue(sin, proben);
            validateObservationCount(arrno, observations,
                                     firstDataLine + proben);

            proben++;
         }
      }
      else
      {
         std::ostringstream s;
         s << "Incorrect expression dimensions at line " << firstDataLine
           << ": header declares " << arrno << " observations, but the first row "
           << "contains " << valueNo << " data fields; expected " << arrno
           << " expression values or " << 2 * arrno
           << " expression/p-value fields.";
         throw s.str();
      }
   }
   catch (std::ios_base::failure& f)
   {
      throw std::string("Could not read data.(Last line empty?)");
   }

   in.exceptions(oldState); //reset excep state
}

//------------------------------------------------------------------------------------

int Microarray_Set::filter(const std::vector<int>& ids, double minMean,
                           double minSigma, int ctlid)
{
   int numIds         = ids.size();
   int numMicroarrays = uarrays.size();
   int numMarkers     = markerset.size();

   for (int i = 0; i < numMarkers; i++)
      markerset[i].isActive = false;

   int numDisabled = numMarkers;

   for (int i = 0; i < numMarkers; i++)
      if (i != ctlid)
      {
         double nx = 0.0, nxx = 0.0;

         for (int j = 0; j < numIds; j++)
         {
            double v = uarrays[ids[j]][i].value;
            nx  += v;
            nxx += v * v;
         }

         double mean     = nx / numMicroarrays;
         double variance = (numMicroarrays * nxx - nx * nx) /
                           (numMicroarrays * numMicroarrays);
         double stdev    = std::sqrt(variance);

         if (mean >= minMean && stdev >= mean * minSigma)
         {
            markerset[i].isActive = true;
            numDisabled--;
         }
      }

   return numDisabled;
}

//------------------------------------------------------------------------------------

int Microarray_Set::filter(double minMean, double minSigma, int ctlid)
{
   int numMicroarrays = uarrays.size();

   std::vector<int> v;

   for (int id = 0; id < numMicroarrays; id++)
      v.push_back(id);

   return filter(v, minMean, minSigma, ctlid);
}

//------------------------------------------------------------------------------------

void Microarray_Set::computeMarkerVariance(const std::vector<int> *arrays)
{
   int numMarkers = markerset.size();

   for (int i = 0; i < numMarkers; i++)
      markerset[i].var = variance(i, arrays);
}

//------------------------------------------------------------------------------------

void Microarray_Set::computeMarkerBandwidth(const std::vector<int> *arrays)
{
   computeMarkerVariance(arrays);

   int n = (arrays ? arrays->size() : uarrays.size());

   double *data = new double[n];

   int numMarkers = markerset.size();

   for (int i = 0; i < numMarkers; i++)
   {
      double prop = 1.06; // Gaussian
      int dim = 1;        // dimension of data

      double stdev = std::sqrt(markerset[i].var);

      for (int j = 0; j < n; j++)
         data[j] = uarrays[arrays ? arrays->at(j) : j][i].value;

      std::sort(data, data + n);

      double iqr = interQuartileRange(data, n);

      double iqrSig = 0.7413 * iqr; // estimate of sigma
      if (iqrSig == 0.0)
         iqrSig = stdev;

      double sig = std::min(stdev, iqrSig);

      markerset[i].bandwidth = prop * sig * std::pow(n, -1.0 / (dim + 4));
   }

   delete[] data;
}

//------------------------------------------------------------------------------------

double Microarray_Set::variance(int m, const std::vector<int> *arrays)
{
   // compute the variance of a gene expression vector

   int n = (arrays ? arrays->size() : uarrays.size());

   double s = 0.0, ss = 0.0;

   for (int i = 0; i < n; i++)
   {
      double v = uarrays[arrays ? arrays->at(i) : i][m].value;
      s  += v;
      ss += v * v;
   }

   return (ss - s * s / n) / (n - 1);
}

//------------------------------------------------------------------------------------

void Microarray_Set::getHighLowPercent(double x, int mId, std::vector<int>& lower,
                                       std::vector<int>& upper)
{
   std::vector<ArrayValuePair> sortArray;
   SortIncreasing_ArrayValuePair sorter;

   int numMicroarrays = uarrays.size();

   for (int id = 0; id < numMicroarrays; id++)
      sortArray.push_back(ArrayValuePair(id, uarrays[id][mId].value));

   std::sort(sortArray.begin(), sortArray.end(), sorter);

   int idPercNo = numMicroarrays * x;

   for (int id = 0; id < idPercNo; id++)
   {
      lower.push_back(sortArray[id].arrayId);
      upper.push_back(sortArray[numMicroarrays - idPercNo - 1 + id].arrayId);
   }
}

//------------------------------------------------------------------------------------

void Microarray_Set::bootStrap(std::vector<int>& boot, const std::vector<int> *arrays)
{
   boot.clear();

   int numIds = (arrays ? arrays->size() : uarrays.size());

   for (int id = 0; id < numIds; id++)
   {
      int r = std::rand() % numIds;
      boot.push_back(arrays ? arrays->at(r) : r);
   }
}

//------------------------------------------------------------------------------------

void Microarray_Set::addNoise()
{
   int numMicroarrays = uarrays.size();
   int numMarkers     = markerset.size();

   for (int id = 0; id < numMicroarrays; id++)
      for (int mid = 0; mid < numMarkers; mid++)
      {
         double r = std::rand();
         double noise = (r / RAND_MAX) * 1e-10;

         uarrays[id][mid].value += noise;
      }
}

//------------------------------------------------------------------------------------

static void BuildRankCache(const Microarray_Set& data, int maNum,
                           const std::vector<int> *arrays,
                           const std::vector<bool>& markersNeeded,
                           std::vector<int>& rankRows,
                           std::vector<int>& rankCache)
{
   // Ranks belong to the exact selected/bootstrap observation sequence.  Position,
   // not original array ID, is the secondary sort key so repeated bootstrap
   // observations retain the same tie behavior as the legacy per-edge sorts.
   if (arrays != NULL && static_cast<int>(arrays->size()) != maNum)
      throw std::string("Rank-cache observation count does not match array selection.");

   const std::size_t numMarkers = data.markerset.size();

   if (markersNeeded.size() != numMarkers)
      throw std::string("Rank-cache marker selection has an invalid size.");

   const std::size_t numRankedMarkers =
      std::count(markersNeeded.begin(), markersNeeded.end(), true);

   if (maNum <= 0 || numRankedMarkers > rankCache.max_size() / maNum)
      throw std::string("Rank-cache dimensions are invalid or too large.");

   rankRows.assign(numMarkers, -1);
   rankCache.resize(numRankedMarkers * static_cast<std::size_t>(maNum));

   if (numRankedMarkers == 0)
      return;

   std::vector<RankedObservation> observations(maNum);
   Sort_RankedObservation sorter;
   int rankRow = 0;

   for (std::size_t marker = 0; marker < numMarkers; marker++)
   {
      if (!markersNeeded[marker])
         continue;

      rankRows[marker] = rankRow++;

      for (int position = 0; position < maNum; position++)
      {
         int arrayId = (arrays ? arrays->at(position) : position);
         observations[position] =
            RankedObservation(data.uarrays[arrayId][marker].value, position);
      }

      std::sort(observations.begin(), observations.end(), sorter);

      std::size_t offset = static_cast<std::size_t>(rankRows[marker]) * maNum;

      for (int rank = 0; rank < maNum; rank++)
         rankCache[offset + observations[rank].position] = rank + 1;
   }
}

//------------------------------------------------------------------------------------

static double Compute_Pairwise_MI(const int *xranks, const int *yranks, int N,
                                  int nparLimit,
                                  AdaptivePartitionWorkspace& workspace)
{
   const int M = nparLimit;

   int npar = 1; maxNpar = 1;
   int run  = 0;

   double xcor = 0.0;

   if (!workspace.matches(N, M))
      throw std::string("Adaptive-partitioning workspace dimensions do not match MI input.");

   workspace.resetEdge();

   std::vector<int>& poc    = workspace.poc;
   std::vector<int>& kon    = workspace.kon;
   std::vector<int>& poradi = workspace.poradi;
   std::vector<int>& marg   = workspace.marg;
   std::vector<int>& apor   = workspace.apor;
   std::vector<unsigned char>& quadrant = workspace.quadrant;
   std::vector<int>& NN     = workspace.NN;
   std::vector<int>& amarg  = workspace.amarg;

   while (npar > 0)
   {
      run++;

      int np   = npar - 1;
      int apoc = poc[np];
      int akon = kon[np];
      int Nex  = akon - apoc + 1;

      for (int i = 0; i < Nex; i++)
         apor[i] = poradi[apoc + i - 1];

      int ave1 = std::floor((marg[np] + marg[np + 2 * M]) / 2);
      int ave2 = std::floor((marg[np + M] + marg[np + 3 * M]) / 2);

      std::fill(NN.begin(), NN.end(), 0);

      for (int i = 0; i < Nex; i++)
      {
         int k = apor[i] - 1;

         int j = (xranks[k] <= ave1 ? 0 : 2) + (yranks[k] <= ave2 ? 0 : 1);

         quadrant[i] = static_cast<unsigned char>(j);
         NN[j]++;
      }

      double c   = Nex / 4.0;
      double sum = 0.0;

      for (int i = 0; i < 4; i++)
      {
         double d = NN[i] - c;
         sum += d * d;
      }

      double tst = 4 * sum / Nex;

      if (tst > 7.8 || run == 1)
      {
         amarg[ 0] = amarg[ 1] = marg[np];
         amarg[ 2] = amarg[ 3] = ave1 + 1;
         amarg[ 4] = amarg[ 6] = marg[np + M];
         amarg[ 5] = amarg[ 7] = ave2 + 1;
         amarg[ 8] = amarg[ 9] = ave1;
         amarg[10] = amarg[11] = marg[np + 2 * M];
         amarg[12] = amarg[14] = ave2;
         amarg[13] = amarg[15] = marg[np + 3 * M];

         // Pack every non-leaf child in one stable pass.  The old code scanned
         // the parent once per child quadrant and copied through a temporary
         // vector; these offsets preserve both quadrant and observation order.
         int writePosition[4] = { -1, -1, -1, -1 };
         int nextPosition = apoc - 1;

         for (int i = 0; i < 4; i++)
            if (NN[i] > 2)
            {
               writePosition[i] = nextPosition;
               nextPosition += NN[i];
            }

         for (int i = 0; i < Nex; i++)
         {
            int childQuadrant = quadrant[i];

            if (writePosition[childQuadrant] >= 0)
               poradi[writePosition[childQuadrant]++] = apor[i];
         }

         npar--;

         for (int i = 0; i < 4; i++)
            if (NN[i] > 2)
            {
               if (++npar > M)
                  throw std::string("Exceeded npar limit!");

               if (npar > maxNpar)
                  maxNpar = npar;

               akon = apoc + NN[i] - 1;

               int np = npar - 1;

               poc[np] = apoc;
               kon[np] = akon;

               for (int j = 0; j < 4; j++)
                  marg[np + j * M] = amarg[i + 4 * j];

               apoc = akon + 1;
            }
            else if (NN[i] > 0)
            {
               double Nx = amarg[i +  8] - amarg[i] + 1;
               double Ny = amarg[i + 12] - amarg[i + 4] + 1;

               xcor += NN[i] * std::log(NN[i] / (Nx * Ny));
            }
      }
      else
      {
         double Nx = marg[np + 2 * M] - marg[np] + 1;
         double Ny = marg[np + 3 * M] - marg[np + M] + 1;

         xcor += Nex * std::log(Nex / (Nx * Ny));

         npar--;
      }
   }

   return (xcor / N + std::log(N));
}

//------------------------------------------------------------------------------------

double Microarray_Set::calculateMI(int maNum, int probeId1, int probeId2,
                                   double threshold, double noise2, int nparLimit,
                                   const std::vector<int>& rankCache,
                                   const std::vector<int>& rankRows,
                                   AdaptivePartitionWorkspace& workspace) const
{
   // compute mutual information between two gene expression vectors;
   // zero is returned if there is no connection

   if (isSameGene(probeId1, probeId2))
      return 0.0;

   if (rankRows[probeId1] < 0 || rankRows[probeId2] < 0)
      throw std::string("Internal error: MI requested for an uncached marker rank.");

   std::size_t offset1 = static_cast<std::size_t>(rankRows[probeId1]) * maNum;
   std::size_t offset2 = static_cast<std::size_t>(rankRows[probeId2]) * maNum;

   double mi = Compute_Pairwise_MI(&rankCache[offset1], &rankCache[offset2],
                                   maNum, nparLimit, workspace);

   if (!std::isfinite(mi))
   {
      std::ostringstream s;
      s << "Non-finite MI calculated for edge "
        << markerset[probeId1].accnum.substr(1) << " -> "
        << markerset[probeId2].accnum.substr(1) << " using " << maNum
        << " observations.";
      throw s.str();
   }

   if (mi < threshold)
      return 0.0;

   if (noise2 <= 0.0)
      return mi;

   double v1 = markerset[probeId1].var;
   double v2 = markerset[probeId2].var;

   double lambda = (v1 / (v1 - noise2)) * (v2 / (v2 - noise2));

   double value  = 1 + (std::exp(2 * mi) - 1) * (1 - 1 / lambda);

   double correctedMi = mi + 0.5 * std::log(value);

   if (!std::isfinite(correctedMi))
   {
      std::ostringstream s;
      s << "Non-finite noise-corrected MI calculated for edge "
        << markerset[probeId1].accnum.substr(1) << " -> "
        << markerset[probeId2].accnum.substr(1) << " using " << maNum
        << " observations.";
      throw s.str();
   }

   return correctedMi;
}

//------------------------------------------------------------------------------------

void Microarray_Set::computeOneRow(int maNum, Matrix& matrix, double threshold,
                                   int row_idx, int numMarkers, int controlId,
                                   bool half_matrix, bool symmetric, double noise2,
                                   int nparLimit,
                                   const std::vector<int>& rankCache,
                                   const std::vector<int>& rankRows,
                                   AdaptivePartitionWorkspace& workspace) const
{
   // this function computes one row of the adjacency matrix; it is called by
   // createEdgeMatrix(); note that since the adjacency matrix is symmetric,
   // only the upper right triangle is computed; so here this function computes only
   // the "upper right triangle" of the row

   for (int j = (half_matrix ? row_idx + 1 : 0); j < numMarkers; j++)
      if (j != controlId && markerset[j].isActive)
      {
         double edge = calculateMI(maNum, row_idx, j, threshold, noise2, nparLimit,
                                   rankCache, rankRows, workspace);
         if (edge != 0.0)
            matrix.addNode(row_idx, j, edge, symmetric);
      }
}

//------------------------------------------------------------------------------------

void Microarray_Set::createEdgeMatrix(int maNum, Matrix& matrix, double threshold,
                                      int controlId, double noise2, int nparLimit,
                                      const std::vector<int>& ids,
                                      const std::vector<int> *arrays) const
{
   // if controlId == -1, there is no constraint (use all arrays to compute the
   // mutual information; if ids.size == 0, all genes will be computed; otherwise,
   // only selected genes will be computed; arrays points to a vector of array ids
   // used for mutual information computation

   std::time_t t1, t2;
   std::time(&t1); // get a timestamp to assess how much time this will take

   int numMarkers = markerset.size();
   int count      = (ids.size() == 0 ? numMarkers : ids.size());
   int step       = std::ceil(0.1 * count);

   std::vector<bool> markersNeeded(numMarkers, false);
   bool hasSource = false;

   if (ids.size() == 0)
   {
      for (int marker = 0; marker < numMarkers; marker++)
         if (marker != controlId && markerset[marker].isActive)
            hasSource = true;
   }
   else
      for (std::vector<int>::const_iterator id = ids.begin(); id != ids.end(); ++id)
         if (*id != controlId)
         {
            markersNeeded[*id] = true; // requested hubs may be inactive
            hasSource = true;
         }

   if (hasSource)
      for (int marker = 0; marker < numMarkers; marker++)
         if (marker != controlId && markerset[marker].isActive)
            markersNeeded[marker] = true;

   std::size_t numRankedMarkers =
      std::count(markersNeeded.begin(), markersNeeded.end(), true);
   double rankCacheMiB = numRankedMarkers * static_cast<double>(maNum) *
                         sizeof(int) / (1024.0 * 1024.0);

   std::vector<int> rankRows;
   std::vector<int> rankCache;
   AdaptivePartitionWorkspace workspace;

   try
   {
      BuildRankCache(*this, maNum, arrays, markersNeeded, rankRows, rankCache);
   }
   catch (const std::bad_alloc&)
   {
      std::ostringstream s;
      s << "Unable to allocate rank cache for " << numRankedMarkers << " genes x "
        << maNum << " observations (" << rankCacheMiB << " MiB).";
      throw s.str();
   }
   catch (const std::length_error&)
   {
      std::ostringstream s;
      s << "Rank cache is too large for " << numRankedMarkers << " genes x "
        << maNum << " observations (" << rankCacheMiB << " MiB).";
      throw s.str();
   }

   std::time(&t2);
   std::cout << "[MI] Cached ranks for " << numRankedMarkers << " genes x " << maNum
             << " observations in " << std::difftime(t2, t1) << " seconds ("
             << rankCacheMiB << " MiB)." << std::endl;

   if (hasSource)
      try
      {
         workspace.initialize(maNum, nparLimit);
      }
      catch (const std::bad_alloc&)
      {
         std::ostringstream s;
         s << "Unable to allocate adaptive-partitioning workspace for " << maNum
           << " observations and partition limit " << nparLimit << ".";
         throw s.str();
      }
      catch (const std::length_error&)
      {
         std::ostringstream s;
         s << "Adaptive-partitioning workspace is too large for " << maNum
           << " observations and partition limit " << nparLimit << ".";
         throw s.str();
      }

   if (ids.size() == 0) // all genes will be computed
   {
      matrix.createEntries(count);

      for (int i = 0; i < count; i++)
      {
         if (i != controlId && markerset[i].isActive)
            computeOneRow(maNum, matrix, threshold, i, numMarkers, controlId,
                          true, true, noise2, nparLimit, rankCache, rankRows,
                          workspace);

         if ((i + 1) % step == 0)
         {
            std::time(&t2);
            std::cout << 10 * (i + 1) / step << "%, time: "
                      << std::difftime(t2, t1) << std::endl;
         }
      }
   }
   else // only selected genes will be computed
      for (int i = 0; i < count; i++)
         if (ids[i] != controlId)
         {
            while (matrix.nmv.size() <= ids[i])
               matrix.nmv.push_back(NodeMap());

            computeOneRow(maNum, matrix, threshold, ids[i], numMarkers, controlId,
                          false, false, noise2, nparLimit, rankCache, rankRows,
                          workspace);

            if ((i + 1) % step == 0)
            {
               std::time(&t2);
               std::cout << 10 * (i + 1) / step << "%, time: "
                         << std::difftime(t2, t1) << std::endl;
            }
         }

   std::time(&t2);
   std::cout << "Gene: " << count << " Time: " << std::difftime(t2, t1) << std::endl;
}
