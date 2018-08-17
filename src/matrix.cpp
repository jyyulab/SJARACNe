//------------------------------------------------------------------------------------
// Copyright (C) 2003  Columbia Genome Center * All Rights Reserved. *
//
// Modifications by S.V. Rice, 2017
//------------------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <functional>
#include <ios>
#include <iterator>
#include <sstream>
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

class GenePair
{
public:
   GenePair()
      : x(), y(), xi(), yi(), maId() { }

   GenePair(double inX, double inY, int inMaId)
      : x(inX), y(inY), xi(), yi(), maId(inMaId) { }

   double x;    // Get_X(), Set_X()
   double y;    // Get_Y(), Set_Y()
   int    xi;   // index of x; Get_XI(), Set_XI()
   int    yi;   // index of y; Get_YI(), Set_YI()
   int    maId; // Get_MaID(), Set_MaID()
};

typedef std::vector<GenePair> GenePairVector;

typedef std::binary_function<GenePair, GenePair, bool> GenePairBinaryFunction;

class Sort_X : GenePairBinaryFunction
{
public:
   bool operator()(const GenePair& a, const GenePair& b) const
   {
      return (a.x < b.x || a.x == b.x && a.maId < b.maId);
   }
};

class Sort_Y : GenePairBinaryFunction
{
public:
   bool operator()(const GenePair& a, const GenePair& b) const
   {
      return (a.y < b.y || a.y == b.y && a.maId < b.maId);
   }
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

void Matrix::read(std::istream& in, Microarray_Set& data, const Parameter& p)
{
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

      while (nmv.size() <= geneId1)
         nmv.push_back(NodeMap());

      std::getline(sin, label, '\t');
      label = "_" + label;

      while (sin.good())
      {
         std::string value;

         std::getline(sin, value, '\t');

         double mi = std::atof(value.c_str());

         if (mi >= p.threshold)
         {
            int geneId2 = data.getProbeId(label);
            if (geneId2 == -1)
               throw "Cannot find marker: " + label + " in the ADJ file!";

            while (nmv.size() <= geneId2)
               nmv.push_back(NodeMap());

            saveNode(geneId1, geneId2, mi);
         }

         std::getline(sin, label, '\t');
         label = "_" + label;
      }

      std::getline(in, line);
   }
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
   int numMarkers = markerset.size();

   for (int i = 0; i < numMarkers; i++)
      if (markerset[i].accnum == accnum)
         return i;

   return -1; // did not find it
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
         readMarkerNoPvalue(fin, proben++);

         while (in.good() && in.peek() != EOF)
         {
            std::getline(in, line);
            std::istringstream sin(line);

            std::getline(sin, accnum, '\t');
            std::getline(sin, label,  '\t');

            accnum = "_" + accnum;

            Marker m(proben, accnum, label);
            Set_Marker(proben, m);

            if (readMarkerNoPvalue(sin, proben) != arrno)
            {
               std::ostringstream s;
               s << "Incorrect data format at line no: " << proben + 2;
               throw s.str();
            }

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
         readMarkerWithPvalue(fin, proben++);

         while (in.good() && in.peek() != EOF)
         {
            std::getline(in, line);
            std::istringstream sin(line);

            std::getline(sin, accnum, '\t');
            std::getline(sin, label,  '\t');

            accnum = "_" + accnum;

            Marker m(proben, accnum, label);
            Set_Marker(proben, m);

            if (readMarkerWithPvalue(sin, proben) != arrno)
            {
               std::ostringstream s;
               s << "Incorrect data format at line no: " << proben + 2;
               throw s.str();
            }

            proben++;
         }
      }
      else
         throw std::string("Incorrect file format: header line doesn't match the "
                           "rest of the data.");
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

static double Compute_Pairwise_MI(GenePairVector& pairs, int nparLimit)
{
   const int M = nparLimit;
   const int N = pairs.size();

   for (int i = 0; i < N; i++)
   {
      pairs[i].xi = i;
      pairs[i].yi = i;
   }

   Sort_X X_Sorter;
   std::sort(pairs.begin(), pairs.end(), X_Sorter);

   std::vector<int> xranks(N);

   for (int i = 0; i < N; i++)
      xranks[pairs[i].xi] = i + 1;

   Sort_Y Y_Sorter;
   std::sort(pairs.begin(), pairs.end(), Y_Sorter);

   std::vector<int> yranks(N);

   for (int i = 0; i < N; i++)
      yranks[pairs[i].yi] = i + 1;

   int npar = 1; maxNpar = 1;
   int run  = 0;

   double xcor = 0.0;

   std::vector<int> poc(M, 1), kon(M, N), poradi(N), marg(4 * M, 0);

   for (int i = 0; i < N; i++)
      poradi[i] = i + 1;

   marg[0] = marg[M] = 1;
   marg[2 * M] = marg[3 * M] = N;

   while (npar > 0)
   {
      run++;

      int np   = npar - 1;
      int apoc = poc[np];
      int akon = kon[np];
      int Nex  = akon - apoc + 1;

      std::vector<int> apor(Nex);

      for (int i = 0; i < Nex; i++)
         apor[i] = poradi[apoc + i - 1];

      int ave1 = std::floor((marg[np] + marg[np + 2 * M]) / 2);
      int ave2 = std::floor((marg[np + M] + marg[np + 3 * M]) / 2);

      std::vector<bool> I(4 * Nex, false);
      std::vector<int>  NN(4, 0);

      for (int i = 0; i < Nex; i++)
      {
         int k = apor[i] - 1;

         int j = (xranks[k] <= ave1 ? 0 : 2) + (yranks[k] <= ave2 ? 0 : 1);

         I[4 * i + j] = true;

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
         std::vector<int> amarg(16);

         amarg[ 0] = amarg[ 1] = marg[np];
         amarg[ 2] = amarg[ 3] = ave1 + 1;
         amarg[ 4] = amarg[ 6] = marg[np + M];
         amarg[ 5] = amarg[ 7] = ave2 + 1;
         amarg[ 8] = amarg[ 9] = ave1;
         amarg[10] = amarg[11] = marg[np + 2 * M];
         amarg[12] = amarg[14] = ave2;
         amarg[13] = amarg[15] = marg[np + 3 * M];

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

               std::vector<int> t;

               for (int j = 0; j < Nex; j++)
                  if (I[i + 4 * j])
                     t.push_back(apor[j]);

               for (int j = apoc - 1; j < akon; j++)
                  poradi[j] = t[j - apoc + 1];

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
                                   const std::vector<int> *arrays) const
{
   // compute mutual information between two gene expression vectors;
   // zero is returned if there is no connection

   if (isSameGene(probeId1, probeId2))
      return 0.0;

   GenePairVector pairs;

   for (int i = 0; i < maNum; i++)
   {
      int index = (arrays ? arrays->at(i) : i);

      double x = uarrays[index][probeId1].value;
      double y = uarrays[index][probeId2].value;

      pairs.push_back(GenePair(x, y, i));
   }

   double mi = Compute_Pairwise_MI(pairs, nparLimit);

   if (mi < threshold)
      return 0.0;

   if (noise2 <= 0.0)
      return mi;

   double v1 = markerset[probeId1].var;
   double v2 = markerset[probeId2].var;

   double lambda = (v1 / (v1 - noise2)) * (v2 / (v2 - noise2));

   double value  = 1 + (std::exp(2 * mi) - 1) * (1 - 1 / lambda);

   return (mi + 0.5 * std::log(value));
}

//------------------------------------------------------------------------------------

void Microarray_Set::computeOneRow(int maNum, Matrix& matrix, double threshold,
                                   int row_idx, int numMarkers, int controlId,
                                   const std::vector<int> *arrays, bool half_matrix,
                                   bool symmetric, double noise2, int nparLimit) const
{
   // this function computes one row of the adjacency matrix; it is called by
   // createEdgeMatrix(); note that since the adjacency matrix is symmetric,
   // only the upper right triangle is computed; so here this function computes only
   // the "upper right triangle" of the row

   for (int j = (half_matrix ? row_idx + 1 : 0); j < numMarkers; j++)
      if (j != controlId && markerset[j].isActive)
      {
         double edge = calculateMI(maNum, row_idx, j, threshold, noise2, nparLimit,
                                   arrays);
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

   if (ids.size() == 0) // all genes will be computed
   {
      matrix.createEntries(count);

      for (int i = 0; i < count; i++)
      {
         if (i != controlId && markerset[i].isActive)
            computeOneRow(maNum, matrix, threshold, i, numMarkers, controlId,
                          arrays, true, true, noise2, nparLimit);

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
                          arrays, false, false, noise2, nparLimit);

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
