//------------------------------------------------------------------------------------
// Copyright (C) 2003, 2008  Columbia Genome Center * All Rights Reserved. *
//
// Modifications by S.V. Rice, 2017
//------------------------------------------------------------------------------------

#ifndef UTIL_H__
#define UTIL_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

//#define M_PI 3.14159265358979323846 (this is defined in cmath)

//Utility Math functions

double median(const std::vector<double>& sortedData, const int n);

double interQuartileRange(const double *sortedData, const int size);

double normalPDF(double dx, double sigma);

double multinormalPDF(double dx, double dy, double sigmaX, double sigmaY);

#endif
