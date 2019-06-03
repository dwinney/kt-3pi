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
#include <cmath>    // complex pow
#include <iostream>

using std::vector;

void gauleg(double x1, double x2, double x[], double w[], int n); // Gaussian-Legendre weights

// Utility functions that take in a vector of complex<doubles> and return vectors of the same size
// containing only the real or imaginary parts.
//-----------------------------------------------------------------------------
std::vector<double> vec_real( std::vector<std::complex<double>> fx);
std::vector<double> vec_imag( std::vector<std::complex<double>> fx);

#endif
