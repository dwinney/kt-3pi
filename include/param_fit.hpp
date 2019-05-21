// Routines to fit a lineshape to a polynomial (z, theta) expansion around center of dalitz plot.
//
// Dependencies: dalitz.cpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _POLY_FIT_
#define _POLY_FIT_

#include "dalitz.hpp"
#include <iomanip>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using std::vector;
using std::setw;

template <class T>
class param_fit : public dalitz<T>
{
protected:
  // Polynomial expansion parameters up to order 5/2 in z
  int n_params = 2;
  double Norm = 1., alpha = 0., beta = 0., gamma = 0., delta = 0.;
  double scale = 1.e-3;

  double F_poly(double s, double t);

  double kin_kernel(double s, double t); // Kinematic Kernel in dalitz region integral
  double chi_squared(const double *par); // Chi-squared between input line-shape and polynomial.

// ---------------------------------------------------------------------------
public:
  param_fit(T my_amp) : dalitz<T>(my_amp){};

  void set_params(int n, const double *par);
  void print_params(int a = 0);

  void extract_params(double N);
// ---------------------------------------------------------------------------
};
// ---------------------------------------------------------------------------
#endif
