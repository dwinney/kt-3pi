// Routines to fit a lineshape to another lineshape by minimizing a chi2-type integral over
// the entire dalitz region
//
// Dependencies: dalitz.hpp, ROOT
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _POLY_FIT_
#define _POLY_FIT_

#include "dalitz.hpp"

#include <string>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using std::string;

// ---------------------------------------------------------------------------
// The dalitz_fit class allows two amplitudes to be fit to each other by minimizing a chi_squared type integral
// over the entire allowed decay region.
//
// my_amp data_amp(&my_inputs);
// dalitz_fit<my_amp, fitting_amp> fit(data_amp);
// fitting_amp finished_fit = fit.extract_params(number_of_variables);
//
// The user must define the function amplitude::set_params(int n, const double *par)
// to properly set the n fitting parameters in the model.
// ---------------------------------------------------------------------------
template <class T, class F>
class dalitz_fit : public dalitz<T>
{
// ---------------------------------------------------------------------------
protected:
  int n_params;
  int nError = 0; // Default print level for info messages in Minuit (0 - 4)

  F * fit_amp;

  double kin_kernel(double s, double t); // Kinematic Kernel in dalitz region integral
  double chi_squared(const double *par); // Chi-squared between input line-shape and polynomial.

// ---------------------------------------------------------------------------
public:
  dalitz_fit(T * my_amp, F * fit_obj)
  : dalitz<T>(my_amp), fit_amp(fit_obj)
  {};

  void print_params(int a = 0);

  void extract_params(int N);

  //Utility to change print level in TMinuit, default is to surpress all messages
  void set_printLevel(int n){nError = n;};
// ---------------------------------------------------------------------------

};
// ---------------------------------------------------------------------------
#endif
