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
#include "poly_exp.hpp"

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
  int nError = 0; // Default print level for info messages in Minuit (0 - 4)

  poly_exp F_poly; // start with norm = 1. and all params = 0.

  double kin_kernel(double s, double t); // Kinematic Kernel in dalitz region integral
  double chi_squared(const double *par); // Chi-squared between input line-shape and polynomial.

// ---------------------------------------------------------------------------
public:
  param_fit(T my_amp) : dalitz<T>(my_amp){
    F_poly.set_decayMass(my_amp.get_decayMass());
  };

  void set_params(int n, const double *par);
  void print_params(int a = 0);

  void extract_params(double N);
  void set_printLevel(int n){nError = n;};
// ---------------------------------------------------------------------------
};
// ---------------------------------------------------------------------------
#endif
