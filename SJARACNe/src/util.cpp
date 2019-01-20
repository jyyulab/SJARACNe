//------------------------------------------------------------------------------------
// Copyright (C) 2003, 2008  Columbia Genome Center * All Rights Reserved. *
//
// Modifications by S.V. Rice, 2017
//------------------------------------------------------------------------------------

#include "util.h"

//------------------------------------------------------------------------------------

double median(const std::vector<double>& sortedData, const int n)
{
   if (n == 0)
      return 0.0;

   const int lhs = (n - 1) / 2;
   const int rhs = n / 2;

   if (lhs == rhs)
      return sortedData[lhs];
   else
      return (sortedData[lhs] + sortedData[rhs]) / 2.0;
}

//------------------------------------------------------------------------------------

double interQuartileRange(const double *sortedData, const int size)
{
   int medianIndex = (size + 2) / 2;

   std::vector<double> subset;

   for (int i = 0; i < medianIndex - 1; i++)
      subset.push_back(sortedData[i]);

   double q1 = median(subset, subset.size());

   medianIndex = (size + 1) / 2;

   subset.clear();

   for (int i = medianIndex; i < size; i++)
      subset.push_back(sortedData[i]);

   double q3 = median(subset, subset.size());

   return q3 - q1;
}

//------------------------------------------------------------------------------------

double normalPDF(double dx, double sigma)
{
  return std::exp(-0.5 * std::pow(dx / sigma, 2)) / (std::sqrt(2 * M_PI) * sigma);
}

//------------------------------------------------------------------------------------

double multinormalPDF(double dx, double dy, double sigmaX, double sigmaY)
{
   return std::exp(-0.5 * (std::pow(dx / sigmaX, 2) + std::pow(dy / sigmaY, 2))) /
          (2 * M_PI * sigmaX * sigmaY);
}
