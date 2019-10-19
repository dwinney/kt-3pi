// Routines to fit a lineshape to another by minimizing a chi2-type integral over
// the entire dalitz region
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
// ---------------------------------------------------------------------------

class dalitz_fit : public dalitz
{
private:
  double smax = dalitz::amp->kinematics.smax();
  double smin = dalitz::amp->kinematics.smin();
  double offset = dalitz::offset;

// ---------------------------------------------------------------------------
protected:
  int n_params;
  int nError = 0; // Default print level for info messages in Minuit (0 - 4)

  amplitude * fit_amp;

  double chi_squared(const double *par); // Chi-squared between input line-shape and polynomial.

// ---------------------------------------------------------------------------
public:
  dalitz_fit(amplitude * my_amp, amplitude * fit_obj)
  : dalitz(my_amp), fit_amp(fit_obj)
  {};

  // print currently stored parameters to the command line
  void print_params(int a = 0);

  // plot the % deviation from resulting fit and the input amplitude
  void plot_deviation();

  // preform fit of N free parameters
  void extract_params(int N);

  //Utility to change print level in TMinuit, default is to surpress all messages
  void set_printLevel(int n){nError = n;};
// ---------------------------------------------------------------------------

};
// ---------------------------------------------------------------------------
#endif
