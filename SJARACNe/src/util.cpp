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

   // Optimization: Use bitwise shift for division by 2 (compiler may optimize, but explicit is clearer)
   // For n > 0, (n-1)/2 and n/2 can be calculated efficiently
   const int lhs = (n - 1) >> 1;  // Equivalent to (n-1)/2, but faster
   const int rhs = n >> 1;        // Equivalent to n/2, but faster

   if (lhs == rhs)
      return sortedData[lhs];
   else
      return (sortedData[lhs] + sortedData[rhs]) * 0.5;  // Multiplication by 0.5 is faster than division by 2.0
}

//------------------------------------------------------------------------------------

double interQuartileRange(const double *sortedData, const int size)
{
   // Optimization: Calculate quartiles directly from sorted data without creating intermediate vectors
   // This eliminates memory allocations and vector operations
   // Matches original logic: Q1 from first (size+2)/2 - 1 elements, Q3 from elements (size+1)/2 onwards
   
   if (size == 0)
      return 0.0;
   
   // Calculate Q1: median of first (medianIndex - 1) elements
   // Original: medianIndex = (size + 2) / 2, then use indices 0 to medianIndex - 2
   int q1End = ((size + 2) / 2) - 1;  // Last index of first half (exclusive in original)
   int q1Size = q1End;  // Number of elements in first half
   
   double q1;
   if (q1Size <= 0)
      q1 = sortedData[0];
   else if (q1Size == 1)
      q1 = sortedData[0];
   else
   {
      // Calculate median of first q1Size elements
      int lhs = (q1Size - 1) >> 1;
      int rhs = q1Size >> 1;
      if (lhs == rhs)
         q1 = sortedData[lhs];
      else
         q1 = (sortedData[lhs] + sortedData[rhs]) * 0.5;
   }
   
   // Calculate Q3: median of elements from (size+1)/2 onwards
   int q3Start = (size + 1) / 2;
   int q3Size = size - q3Start;
   
   double q3;
   if (q3Size <= 0)
      q3 = sortedData[size - 1];
   else if (q3Size == 1)
      q3 = sortedData[q3Start];
   else
   {
      // Calculate median of last q3Size elements
      int lhs = (q3Size - 1) >> 1;
      int rhs = q3Size >> 1;
      if (lhs == rhs)
         q3 = sortedData[q3Start + lhs];
      else
         q3 = (sortedData[q3Start + lhs] + sortedData[q3Start + rhs]) * 0.5;
   }
   
   return q3 - q1;
}

//------------------------------------------------------------------------------------

double normalPDF(double dx, double sigma)
{
   // Optimization: Replace std::pow(x, 2) with x * x (faster)
   // Pre-compute constant: 1 / sqrt(2 * PI) for better performance
   static const double INV_SQRT_2PI = 1.0 / std::sqrt(2.0 * M_PI);
   
   double dx_sigma = dx / sigma;
   return std::exp(-0.5 * dx_sigma * dx_sigma) * INV_SQRT_2PI / sigma;
}

//------------------------------------------------------------------------------------

double multinormalPDF(double dx, double dy, double sigmaX, double sigmaY)
{
   // Optimization: Replace std::pow(x, 2) with x * x (faster)
   // Pre-compute constant: 1 / (2 * PI) for better performance
   static const double INV_2PI = 1.0 / (2.0 * M_PI);
   
   double dx_sigmaX = dx / sigmaX;
   double dy_sigmaY = dy / sigmaY;
   return std::exp(-0.5 * (dx_sigmaX * dx_sigmaX + dy_sigmaY * dy_sigmaY)) *
          INV_2PI / (sigmaX * sigmaY);
}
