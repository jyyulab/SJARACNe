//------------------------------------------------------------------------------------
// Copyright (C) 2003, 2004  Columbia Genome Center * All Rights Reserved. *
//
// Modifications by S.V. Rice, 2017
//------------------------------------------------------------------------------------

#include <cctype>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include "param.h"

//------------------------------------------------------------------------------------

const double Parameter::default_threshold  =  0.00; // mi threshold
const double Parameter::default_pvalue     =  1.00; // p-value for mi threshold
const double Parameter::default_eps        =  1.00; // DPI tolerance
const double Parameter::default_sigma      = 99.00; // kernel width
const int    Parameter::default_sample     =  0;    // sample number
const double Parameter::default_percent    =  0.35; // high/low percentage
const double Parameter::default_mean       =  0.00; // filter mean
const double Parameter::default_cv         =  0.00; // filter standard deviation
const double Parameter::default_correction =  0.00; // array measurement noise level
const int    Parameter::default_nparLimit  = 20;    // max allowed value of npar
const int    Parameter::default_seed       = 1;     // Initial seed for random number generator
//------------------------------------------------------------------------------------

bool equalIgnoreCase(std::string a, std::string b)
{
   int len = a.length();
   if (b.length() != len)
      return false;

   int i = 0;
   while (i < len && std::toupper(a[i]) == std::toupper(b[i]))
      i++;

   return (i == len);
}

//------------------------------------------------------------------------------------

void checkParameter(Parameter &p)
{
   if (p.infile == "")
      throw std::string("No input file specified!");

   if (p.hub != "")
      p.hub = "_" + p.hub;

   if (p.subnetfile != "" && p.hub != "")
      throw std::string("Either supply one hub gene by '-h' or multiple genes in a "
                        "file by '-s', but not both!");

   if (p.condition != "+" && p.condition != "-" && p.condition != "")
      throw std::string("Condition must be '+' or '-'!");

   if ((p.condition == "+" || p.condition == "-") && p.controlId == "")
      throw std::string("Control gene ID must be specified using '-c'!");

   if (p.sigma != Parameter::default_sigma && (p.sigma <= 0.0 || p.sigma >= 1.0))
      throw std::string("Kernel width '-k' must be within (0,1)!");

   if (p.threshold < 0.0)
      throw std::string("MI threshold '-t' must be nonnegative!");

   if (p.threshold > 0.0 && p.pvalue != 1.0)
      std::cout << "P-value will not be used, since a threshold has been specified."
                << std::endl;

   if (p.pvalue <= 0.0 || p.pvalue > 1.0)
      throw std::string("P-value '-p' must be in the range (0,1]!");

   if (p.eps < 0.0 || p.eps > 1.0)
      throw std::string("DPI tolerance '-e' must be within [0,1]!");

   if (p.percent <= 0.0 || p.percent >= 1.0)
      throw std::string("Percentage microarray must be within (0,1)!");

   if (p.mean < 0.0)
      throw std::string("Gene filter mean must be nonnegative!");

   if (p.cv < 0.0)
      throw std::string("Gene filter cv (coefficient of variance) must be "
                        "nonnegative!");

   if (p.correction < 0.0)
      throw std::string("Array measurement noise level '-n' must be nonnegative!");

   if (!equalIgnoreCase(p.verbose, "on") && !equalIgnoreCase(p.verbose, "off"))
      throw std::string("Verbose '-v' must be 'on' or 'off'!");

   if (p.nparLimit < 1)
      throw std::string("Maximum allowed value of npar must be positive!");

   if (p.home_dir != "./")
   {
      int len = p.home_dir.length();
      int b   = p.home_dir.find_last_of("/");

      if (b == std::string::npos || b < len - 1)
         p.home_dir += "/";
   }
}

//------------------------------------------------------------------------------------
// readProbeList() reads the list of nodes to be included in constructing a subnetwork

static int readProbeList(const std::string& infilename,
                         std::vector<std::string>& probe_list)
{
   std::ifstream in(infilename.c_str());
   if (!in.is_open())
      throw "Unable to open " + infilename;

   std::string line;
   int lnum = 0;

   while (in.good() && in.peek() != EOF && in.peek() != '\012')
   {
      std::getline(in, line);
      std::istringstream sin(line);
      std::string gid;
      std::getline(sin, gid);
      gid = "_" + gid;
      probe_list.push_back(gid);
      lnum++;
   }

   in.close();

   return lnum;
}

//------------------------------------------------------------------------------------

void displayParameter(Parameter &p)
{
   std::cout << std::endl;

   std::cout << "[PARA] Input file:    " << p.infile  << std::endl;
   std::cout << "[PARA] Output file:   " << p.outfile << std::endl;

   if (p.threshold > 0.0)
      std::cout << "[PARA] MI threshold:  " << p.threshold << std::endl;
   else
      std::cout << "[PARA] MI P-value:    " << p.pvalue    << std::endl;

   std::cout << "[PARA] DPI tolerance: " << p.eps << std::endl;

   if (p.correction > 0.0)
      std::cout << "[PARA] Correction for MI estimation (array noise level: "
                << p.correction << ")" << std::endl;

   if (p.subnetfile != "")
      std::cout << "[PARA] Subset of probes to reconstruct: " << p.subnetfile
                << " (" << readProbeList(p.subnetfile, p.subnet) << ")" << std::endl;

   if (p.hub != "")
      std::cout << "[PARA] Hub probe to reconstruct: " << p.hub << std::endl;

   if (p.controlId != "")
   {
      std::cout << "[PARA] Control gene:  " << p.controlId << std::endl;
      std::cout << "[PARA] Condition:     " << p.condition << std::endl;
      std::cout << "[PARA] Percentage:    " << p.percent   << std::endl;
   }

   if (p.annotfile != "")
      std::cout << "[PARA] TF annotation list: " << p.annotfile
                << " (" << readProbeList(p.annotfile, p.tf_list) << ")" << std::endl;

   if (p.mean != 0.0 || p.cv != 0.0)
   {
      std::cout << "[PARA] Filter mean:   " << p.mean << std::endl;
      std::cout << "[PARA] Filter CV:     " << p.cv   << std::endl;
   }

   std::cout << "[PARA] Npar limit:    " << p.nparLimit << std::endl;
}

//------------------------------------------------------------------------------------

static std::string getFileName(const std::string& matrixName)
{
   int b = matrixName.find_last_of("/");
   if (b == std::string::npos)
      b = matrixName.find_last_of("\\");

   std::string dirname(""), basename(matrixName);

   if (b != std::string::npos)
   {
      // Extract the directory and the filename if path is included
      basename = matrixName.substr(b + 1);
      dirname  = matrixName.substr(0, b) + "/";
   }

   int c = basename.find_last_of(".");
   if (c != std::string::npos)
       basename = basename.substr(0, c);

   return dirname + basename;
}

//------------------------------------------------------------------------------------

void createOutfileName(Parameter &p)
{
   std::string filename = getFileName(p.infile);

   if (p.hub != "")
      filename += "_h" + p.hub;

   if (p.controlId != "")
      filename += "_c" + p.controlId + (p.condition == "+" ? "H" : "L");

   char buffer[20];

   std::sprintf(buffer, "%0.3g", p.sigma);
   filename += std::string("_k") + buffer;

   if (p.threshold > 0.0)
   {
      std::sprintf(buffer, "%0.2g", p.threshold);
      filename += std::string("_t") + buffer;
   }

   if (p.eps < 1.0)
   {
      std::sprintf(buffer, "%0.2g", p.eps);
      filename += std::string("_e") + buffer;
   }

   if (p.sample > 0)
   {
      std::sprintf(buffer, "%03i", p.sample);
      filename += std::string("_r") + buffer;
   }

   p.outfile = filename + ".adj";
}
