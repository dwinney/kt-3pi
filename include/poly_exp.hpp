// Simple class with a polynomial amplitude mimicing the expansion around center of dalitz plot to test fitting and parameter extraction.
//
// Dependencies: decay_kinematics
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _POLY_
#define _POLY_

#include "decay_kinematics.hpp"
#include "amplitude.hpp"
#include "dalitz.hpp"

#include <iomanip>

using std::setw;

class poly_exp : public amplitude
{
protected:
  double scale = 1.e-3;

  double alpha = 0., beta = 0., gamma = 0., delta = 0.;
  double Norm = 1;

  double dNorm = 0.;
  double dalpha = 0., dbeta = 0., dgamma = 0., ddelta = 0.;

public:
  // Default constructor
  poly_exp(decay_kinematics dec)
  : amplitude(dec)
  {};

  // Constructor with number of params and an array
  poly_exp(int n, const double * par, decay_kinematics dec)
  : amplitude(dec)
  {
    set_params(n, par);
  };

  // Evaluate amplitude
  complex<double> eval(double s, double t);

  bool ERROR = false;
  double error_func(double s, double t);
  void normalize(double gamma_exp);

  // Set and Print
  void set_params(int n, const double * par);
  void set_errors(int n, const double * par);
  void print_params(int a = 0);
};

#endif
