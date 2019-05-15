// General purpose classes for Dalitz plot generation.
// for fitting to polynomial expansion see poly_param_fit.hpp
//
// Dependencies: amp.hpp, aux_math.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _DALITZ_
#define _DALITZ_

#include "constants.hpp"
#include <cmath>
#include<complex>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

using std::complex;
using std::setw;
using std::cout;
using std::endl;

//-----------------------------------------------------------------------------
//  Define a Dalitz plot object with a pointer to the amplitude.
//  The amp should be a complex<double> and a function of 2 doubles, Mandelstam s and t.
//
//  complex<double> (*ptr_amp) (double, double);
//  ptr_func = my_amp;
//
//  dalitz my_dalitz(ptr_amp);
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
  dalitz()
  {
    cout << " Need pointer to amplitude function in dalitz object declaration. Quitting... \n";
    std::exit(1);
  };
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

  // Double differential cross section
  double d2Gamma(double s, double t);

  void plot(const char *n = "");
};

//-----------------------------------------------------------------------------
#endif
