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

#include <iomanip>

using std::setw;

class poly_exp
{
protected:
  double alpha = 0., beta = 0., gamma = 0., delta = 0.;
  double Norm = 1.;
  double scale = 1.e-3;

public:
  // Default constructor
  poly_exp(decay_kinematics dec)
  : kinematics(dec)
  {};

  // Constructor with number of params and an array
  poly_exp(int n, const double * par, decay_kinematics dec)
  : kinematics(dec)
  {
    set_params(n, par);
  };

  decay_kinematics kinematics;

  // Evaluate amplitude squared
  complex<double> eval(double s, double t);

  // Set and Print
  void set_params(int n, const double *par);
  void print_params(int a = 0);
};

#endif
