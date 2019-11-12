// Misc utility functions and classes
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#ifndef _UTIL_
#define _UTIL_

#include <complex>
#include <vector>
#include <iterator>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include <Math/Interpolator.h>
#include <TH2.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TError.h>

using std::string;
using std::vector;
using std::cout;
using std::complex;

//-----------------------------------------------------------------------------
// given int n, outputs right english string
// "1st" vs "2nd" vs "3rd" etc...
// :)
std::string st_nd_rd(int n);

//-----------------------------------------------------------------------------
// Numerical integration routine. Generates Gauss-Legendre weights and abscissas
void gauleg(double x1, double x2, double x[], double w[], int n);

//-----------------------------------------------------------------------------
// Utility functions that take in a vector of complex<doubles> and return vectors of the same size
// containing only the real or imaginary parts.
std::vector<double> vec_real( std::vector<std::complex<double>> fx);
std::vector<double> vec_imag( std::vector<std::complex<double>> fx);

//-----------------------------------------------------------------------------
// Simple function to call ROOT to print a plot of real and imag parts
void quick_plot(vector<double> s, vector<complex<double>> fx, string filename);

// make a symmetric dalitz plot from a file
// optional parameter is whether to plot using white color scheme to see % difference
void quick_dplot(string file, bool DEV = false);

// Print a complex function to a .dat file
void quick_print(vector<double> s, vector<complex<double>> fx, string file);

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
