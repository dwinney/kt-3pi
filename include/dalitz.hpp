// General purpose classes for Dalitz plot generation.
// for fitting to polynomial expansion see poly_param_fit.hpp
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
double offset = 0.00001;

// Center of the Dalitz plot in s and t
double s_c()
  {
    return  (amp.mDec * amp.mDec + 3. * mPi* mPi) / 3.;
  };
double t_c()
{
    return (3.*mPi*mPi + amp.mDec* amp.mDec - s_c())/2.;
};

//-----------------------------------------------------------------------------
public:
  // Default Constructor
  dalitz(T& my_amp) : amp(my_amp){};
//--------------------------------------------------------------------------
  // Lorentz Invariant dimensionless parameters
  double d_norm(){
    return 1. / (2. * amp.mDec * (amp.mDec - 3. * mPi));
  }
  double x(double s, double t);
  double y(double s, double t);

  double x_polar(double z, double theta);
  double y_polar(double z, double theta);

  double z(double s, double t);
  double theta(double s, double t);
//--------------------------------------------------------------------------
  // Double differential cross section
  double d2Gamma(double s, double t);
//--------------------------------------------------------------------------
  void plot(const char *n = "");
};

//-----------------------------------------------------------------------------
#endif
