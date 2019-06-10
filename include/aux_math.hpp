// Auxilary Math Equations
//
// Dependencies: None
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#ifndef _MATH_
#define _MATH_

#include <complex>
#include <vector>
#include <iterator>
#include <cmath>
#include <iostream>

#include "Math/Interpolator.h"

using std::vector;
using std::complex;

void gauleg(double x1, double x2, double x[], double w[], int n); // Gaussian-Legendre weights

//-----------------------------------------------------------------------------
// Utility functions that take in a vector of complex<doubles> and return vectors of the same size
// containing only the real or imaginary parts.
std::vector<double> vec_real( std::vector<std::complex<double>> fx);
std::vector<double> vec_imag( std::vector<std::complex<double>> fx);

//-----------------------------------------------------------------------------
// Wrapper class to better interface with ROOT's interpolation class
// for both real and imaginary parts
class interpolation
{
protected:
  vector<double> s, r_fx, i_fx;
  ROOT::Math::Interpolator r_inter, i_inter;

public:
  // The number of interpolation points used  can be changed here
  const static int N_interp = 200;
  interpolation(){};
  
  interpolation(vector<double> x, vector<complex<double>> fx)
  : s(x), r_fx(vec_real(fx)), i_fx(vec_imag(fx)),
    r_inter(x, vec_real(fx), ROOT::Math::Interpolation::kCSPLINE),
    i_inter(x, vec_imag(fx), ROOT::Math::Interpolation::kCSPLINE)
  { };

  // Copy Constructor (needed for vector<iterations>)
  // Additionally, ROOT's Interpolator object cannot be copied so this will create new
  // Interpolator with the same data.
  interpolation(const interpolation &old_interp)
  : s(old_interp.s), r_fx(old_interp.r_fx), i_fx(old_interp.i_fx),
    r_inter(old_interp.s, old_interp.r_fx, ROOT::Math::Interpolation::kCSPLINE),
    i_inter(old_interp.s, old_interp.i_fx, ROOT::Math::Interpolation::kCSPLINE)
    {};

  // The output, at given s, outputs the real and imaginary parts from the interpolations
  complex<double> operator ()(double s);
};

#endif
