// Classes for Dalitz plot generation and fitting to polynomial expansion with dalitz plot parameters.
//
//
// Dependencies: amp.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _DALITZ_
#define _DALITZ_

#include "amp.hpp"
#include <iostream>
#include <complex>
#include <cmath>
using std::complex;
using std::cout;

//-----------------------------------------------------------------------------
//  Define a Dalitz plot object with a pointer to the amplitude.
//  The amp should be a complex<double> and a function of 2 doubles, Mandelstam s and t.
//
//  complex<double> (*ptr_amp) (double, double);
//  ptr_func = my_amp;
//
//  dalitz my_dalitz(ptr_amp);
//-----------------------------------------------------------------------------

class dalitz : amplitude
{
protected:
complex<double> (*amp) (double, double);
double s_c = (Mass * Mass + 3. * mPi*mPi) / 3.;
double temp_num, temp_den, temp_pref, amp_sqr;

// Polynomial expansion parameters
double Norm, alpha, beta, gamma, delta;

//-----------------------------------------------------------------------------
public:
  dalitz()
  {
    cout << " Need pointer to amplitude function in dalitz object declaration. Quitting... \n";
    std::exit(1);
  };
  dalitz(complex<double> (*my_amp) (double, double))
  {
    amp = my_amp;
  };
//-----------------------------------------------------------------------------
  // Lorentz Invariant dimensionless parameters
  double x(double s, double t);
  double y(double s, double t);

  double x_polar(double z, double theta);
  double y_polar(double z, double theta);

  // Double differential cross section
  double d2Gamma(double s, double t);

  // Polynomial expansion around the center of dalitz plot
  complex<double> F_poly(double z, double theta);
//-----------------------------------------------------------------------------
};

//-----------------------------------------------------------------------------
#endif
