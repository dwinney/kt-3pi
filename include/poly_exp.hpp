// Simple class with a polynomial amplitude mimicing the expansion around center of dalitz plot to test fitting and parameter extraction.
//
// Dependencies: amp.cpp
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

class poly_exp : public decay_kinematics
{
protected:
  double alpha = 0., beta = 0., gamma = 0., delta = 0.;
  double Norm = 1.;
  double scale = 1.e-3;

public:
  // Default constructor
  poly_exp(){};

  // Constructor specifying parameters
  poly_exp(double norm, double alp, double bet, double gam, double del)
          : Norm(norm), alpha(alp), beta(bet), gamma(gam), delta(del) {};

  // Constructor with number of params and an array
  poly_exp(int n, const double * par)
  {
    set_params(n, par);
  }
  // Evaluate amplitude squared
  complex<double> operator ()(double s, double t);

  // Set and Print
  void set_params(int n, const double *par);
  void print_params(int a = 0);
};

#endif
