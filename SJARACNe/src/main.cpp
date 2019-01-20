//------------------------------------------------------------------------------------
// Copyright (C) 2003  Columbia Genome Center * All Rights Reserved. *
//
// Modifications by S.V. Rice, 2017; this version assumes adaptive_partitioning
//------------------------------------------------------------------------------------

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include "matrix.h"
#include "parseargs.h"

//------------------------------------------------------------------------------------

const int NUM_OPTIONS = 20;

const char *option[NUM_OPTIONS] =
{
"-a <algorithm>     default: adaptive_partitioning",
"-c <+/-probeId %>  Conditional network reconstruction, default: NONE [****]\n"
"                   [format: \"+24 0.35\", \"-1973_s_at 0.4\"]",
"-e <tolerance>     DPI tolerance, default: 1",
"-f <mean> <cv>     Gene filter by the mean and coefficient of variance (cv) of\n"
"                   the expression values, default: mean=0, cv=0",
"-H <ARACNE_HOME>   Directory containing ARACNE configuration files,\n"
"                   default: current working directory",
"-h <probeId>       Hub gene (only MI w/ hub gene will be computed),\n"
"                   default: NONE",
"--h | --help       Display this help and exit",
"-i <file>          Input gene expression profile dataset (required)",
"-j <file>          Existing adjacency matrix (.adj) file",
"-k <kernel_width>  Gaussian kernel width (accurate method only),\n"
"                   default: determined by program",
"-l <file>          File containing a list of probes annotated as transcription\n"
"                   factors in the input dataset, default: NONE [***]",
"-N <npar_limit>    Maximum allowed value of npar, default: 20",
"-S <Seed>	    Initial seed for random number generator, default: 1",
"-n <level>         Array measurement noise level, default: 0",
"-o <file>          Output file name (optional) [*]",
"-p <p-value>       P-value for MI threshold (e.g., 1e-7), default: 1 [**]",
"-r <sample_number> Bootstrap sample number, default: 0",
"-s <file>          File containing a list of probes for which a subnetwork will\n"
"                   be constructed, default: NONE",
"-t <threshold>     MI threshold, default: 0",
"-v <verbose>       on|off, default: off"
};

const int NUM_USAGE_NOTES = 4;

const char *usageNotes[NUM_USAGE_NOTES] =
{
"   [*] If no output file is specified by the user, an output will be\n"
"       automatically generated in the same directory as the input file by\n"
"       appending some of the parameter values, such as kernel width, MI\n"
"       threshold, tolerance and so on, at the end of the input file name, and\n"
"       changing the file extension to \".adj\".",
"  [**] If the '-t' option is supplied, it will enforce the program to use the\n"
"       specified MI threshold, therefore the '-p' option will be ignored.\n"
"       Otherwise, the program will automatically determine the MI threshold\n"
"       given the p-value. The default, p-value=1, will preserve all pairwise MI.",
" [***] This option is ideal for transcriptional network reconstruction. If\n"
"       provided, DPI will not remove any connection of a transcription factor (TF)\n"
"       by connections between two probes not annotated as TFs. This option is\n"
"       often used in conjunction with '-s', which specifies a list of probes that\n"
"       are either the same or a subset of the probes specified by '-l'.",
"[****] Conditional network reconstructs the network given a specified probe being\n"
"       most expressed or least expressed. In the format that follows '-c', probeId\n"
"       indicates the probe to be conditioned on; '+' or '-' specify whether the\n"
"       upper or lower tail of the probe's expression should be used as the\n"
"       condition, and '%' is a percentage between (0,1) specifying the proportion\n"
"       of samples used as the conditioning subset. Example usage: \"-c +24 0.35\",\n"
"       \"-c -1973_s_at 0.4\"."
};

//------------------------------------------------------------------------------------
// writeUsage() writes usage to stdout

void writeUsage(const char *progname)
{
   std::printf("Usage: %s [OPTION] ...\n\n", progname);

   for (int i = 0; i < NUM_OPTIONS; i++)
      std::printf("    %s\n", option[i]);

   std::printf("\n");

   for (int i = 0; i < NUM_USAGE_NOTES; i++)
      std::printf("%s\n", usageNotes[i]);
}

//------------------------------------------------------------------------------------

Parameter parseParameter(int argc, char *argv[])
{
   Parameter p; // initializes parameters to default values

   std::string temp;

   // parse arguments using macros defined in parseargs.h
   ARGBEGIN
   {
      case 'a': temp         = ARGF(); break;            // algorithm
      case 'c': temp         = ARGF();                   // condition
                p.condition  = temp.substr(0,1);
                p.controlId  = temp.substr(1);
                p.percent    = std::atof(ARGF());
            break;
      case 'e': p.eps        = std::atof(ARGF()); break; // DPI tolerance
      case 'f': p.mean       = std::atof(ARGF());        // mean
                p.cv         = std::atof(ARGF());        // coefficient of variance
                break;
      case 'H': p.home_dir   = ARGF(); break;            // ARACNE_HOME
      case 'h': p.hub        = ARGF(); break;            // hub gene
      case 'i': p.infile     = ARGF(); break;            // input file
      case 'j': p.adjfile    = ARGF(); break;            // adjacency matrix file
      case 'k': p.sigma      = std::atof(ARGF()); break; // gaussian kernel width
      case 'l': p.annotfile  = ARGF(); break;            // TF annotation file
      case 'N': p.nparLimit  = std::atoi(ARGF()); break; // max npar value
      case 'S': p.seed       = std::atoi(ARGF()); break; // seed
      case 'n': p.correction = std::atof(ARGF()); break; // correction for noise
      case 'o': p.outfile    = ARGF(); break;            // output file
      case 'p': p.pvalue     = std::atof(ARGF()); break; // p-value
      case 'r': p.sample     = std::atoi(ARGF()); break; // bootstrap sample number
      case 's': p.subnetfile = ARGF(); break;            // subset of probes
      case 't': p.threshold  = std::atof(ARGF()); break; // mi threshold
      case 'v': p.verbose    = ARGF(); break;            // verbose
      default : throw std::string("unknown parameter ") + ARGC();
   }
   ARGEND;

   checkParameter(p);

   std::cout << "Displaying parameters" << std::endl;
   displayParameter(p);

   return p;
}

//------------------------------------------------------------------------------------

void findThreshold(int n, Parameter& p)
{
   std::string filename = p.home_dir + "config_threshold.txt";

   std::ifstream infile(filename.c_str());
   if (!infile.is_open())
      throw "Unable to open " + filename;

   std::string line;

   std::getline(infile, line);
   while (line.length() > 0 && line[0] == '>')
      std::getline(infile, line);

   std::istringstream sin(line);

   if (!sin.good() || sin.peek() == EOF)
      throw "Configuration file format error: " + filename;

   double alpha, beta, gamma;
   sin >> alpha >> beta >> gamma;
   p.threshold = (alpha - std::log(p.pvalue)) / (-beta - gamma * n);

   infile.close();
}

//------------------------------------------------------------------------------------

void runStandard(int argc, char *argv[])
{
   Parameter p = parseParameter(argc, argv);
   std::srand(p.seed);
   Microarray_Set data;

   data.read(p.infile);

   if (p.mean > 0.0 || p.cv > 0.0)
      std::cout << data.filter(p.mean, p.cv)
                << " markers disabled due to lack of dynamic range."
                << std::endl << std::endl;

   std::vector<int> lower, upper, *arrays = NULL;

   int controlId = -1;
   int nsample   = data.uarrays.size();

   if (p.controlId != "")
   {
      controlId = data.getProbeId(p.controlId);
      if (controlId == -1)
         throw "Cannot find marker: " + p.controlId;

      data.getHighLowPercent(p.percent, controlId, lower, upper);
      arrays = (p.condition == "+" ? &upper : &lower);
      nsample = arrays->size();
   }

   std::cout << "Marker No: " << data.markerset.size()
             << " (" << data.Get_Num_Active_Markers() << " active)"
             << ", Array No: " << nsample << std::endl;

   if (p.threshold == 0.0 && p.pvalue != 1.0)
   {
      findThreshold(nsample, p);

      std::cout << "MI threshold determined for p=" << p.pvalue << ": " << p.threshold
                << std::endl;
   }

   if (p.correction != 0.0)
      data.computeMarkerVariance(arrays);

   data.computeMarkerBandwidth(arrays);

   std::vector<int> ids;

   if (p.hub != "")
   {
      std::string hub = p.hub;

      if (hub.length() == 0 || hub[0] != '_')
         hub = "_" + hub;

      int hubId = data.getProbeId(hub);
      if (hubId == -1)
         throw "Cannot find the hub probe: " + p.hub + ", nothing to be computed!";

      ids.push_back(hubId);
   }

   int numSubnets = p.subnet.size();

   for (int i = 0; i < numSubnets; i++)
   {
      int gid = data.getProbeId(p.subnet[i]);
      if (gid == -1)
         std::cout << "Cannot find probe: " << p.subnet[i] << " in \""
                   << p.subnetfile << "\" ... ignored." << std::endl;
      else
         ids.push_back(gid);
   }

   Transfac transfac;

   int numTFs = p.tf_list.size();

   for (int i = 0; i < numTFs; i++)
   {
      int gid = data.getProbeId(p.tf_list[i]);
      if (gid == -1)
         std::cout << "Cannot find probe: " << p.tf_list[i] << " in \""
                   << p.annotfile << "\" ... ignored." << std::endl;
      else
     transfac[gid] = 1;
   }

   Matrix matrix;

   if (p.adjfile != "")
      matrix.read(data, p);
   else
   {
      std::vector<int> bs;

      if (p.sample > 0)
      {
         data.bootStrap(bs, arrays);
     arrays = &bs;
      }

      data.addNoise();

      data.createEdgeMatrix(nsample, matrix, p.threshold, controlId, p.correction,
                            p.nparLimit, ids, arrays);
   }

   if (p.eps != 1.0)
   {
      std::cout << "[NETWORK] Applying DPI ..." << std::endl;
      matrix.reduce(p.eps, ids, transfac);
   }

   if (p.outfile == "")
      createOutfileName(p);

   matrix.write(data, ids, p);
}

//------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
   bool displayUsage = false;

   if (argc < 2)
      displayUsage = true;
   else
   {
      std::string firstOption = argv[1];

      if (firstOption == "--help" || firstOption == "--h")
         displayUsage = true;
   }

   if (displayUsage)
   {
      writeUsage(argv[0]);
      return 1;
   }

   try
   {
      runStandard(argc, argv);
   }
   catch (const std::string& s)
   {
      std::cerr << argv[0] << ": " << s << std::endl;
      return 1;
   }

   return 0;
}
