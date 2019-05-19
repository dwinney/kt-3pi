// General purpose classes for Dalitz plot generation.
// for fitting to polynomial expansion see param_fit.hpp
//
// Dependencies: constant.cpp aux_math.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _DALITZ_
#define _DALITZ_

#include "constants.hpp"
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "TGraph2D.h"
#include "TCanvas.h"
#include "TStyle.h"

using std::complex;
using std::setw;
using std::cout;
using std::endl;

//-----------------------------------------------------------------------------
// Define a Dalitz plot object with some other amplitude object containing a model.
// Object argument must inherit from 'amplitude' class and be callable with the signature:
// with the signature:
//
//  complex<double> operator() (double, double)
//
// Initiate a dalitz plot object with:
// model_amp my_model;
// dalitz<model_amp> my_dalitz(my_model);
//-----------------------------------------------------------------------------

template <class T>
class dalitz
{
protected:
T amp;
double normalization = 32. * pow(2.* M_PI * amp.mDec, 3.);
double offset = 0.00001;

//-----------------------------------------------------------------------------
public:
  // Default Constructor
  dalitz(T& my_amp) : amp(my_amp){};

//--------------------------------------------------------------------------
  // Double differential cross section
  double d2Gamma(double s, double t);

//--------------------------------------------------------------------------
  void plot(const char *n = "");
};

//-----------------------------------------------------------------------------
#endif
