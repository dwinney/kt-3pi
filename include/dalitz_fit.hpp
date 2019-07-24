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

using std::string;

// ---------------------------------------------------------------------------
// The dalitz_fit class allows two amplitudes to be fit to each other by minimizing a chi_squared type integral
// over the entire allowed decay region.
//
// The amplitudes for both data and fitting must contain a decay_kinematics object as a
// publically accessable data member.
//
// The user must define the function
// void set_params(int n, const double *par) [to properly set the n fitting parameters in the model]
// and
// complex<double> eval(doube s, double t) [evaluate amplitude at mandelstam s and t]
// ---------------------------------------------------------------------------
template <class T, class F>
class dalitz_fit : public dalitz<T>
{
private:
  double smax = dalitz<T>::amp->kinematics.smax();
  double smin = dalitz<T>::amp->kinematics.smin();
  double offset = dalitz<T>::offset;

// ---------------------------------------------------------------------------
protected:
  int n_params;
  int nError = 0; // Default print level for info messages in Minuit (0 - 4)

  F * fit_amp;

  double chi_squared(const double *par); // Chi-squared between input line-shape and polynomial.

// ---------------------------------------------------------------------------
public:
  dalitz_fit(T * my_amp, F * fit_obj)
  : dalitz<T>(my_amp), fit_amp(fit_obj)
  {};

  void print_params(int a = 0);
  void print_deviation(string filename = "");

  void extract_params(int N);

  //Utility to change print level in TMinuit, default is to surpress all messages
  void set_printLevel(int n){nError = n;};
// ---------------------------------------------------------------------------

};
// ---------------------------------------------------------------------------
#endif
