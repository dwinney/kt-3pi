// Simple class with a polynomial amplitude mimicing the expansion around center of dalitz plot to test fitting and parameter extraction.
//
// Dependencies: amp.cpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _TEST_
#define _TEST_

#include "amplitude.hpp"

class test_poly : public amplitude
{

protected:
  double alpha, beta, gamma, delta;
  double Norm = 1.;
  double scale = 1.e3;

public:
  // Default constructor only two parameters for now
  test_poly(double norm, double alp = 0, double bet = 0, double gam = 0, double del = 0)
          : Norm(norm), alpha(alp), beta(bet), gamma(gam), delta(del) {};

  complex<double> operator ()(double s, double t)
    {
      double zs = z(s,t);
      double thetas = theta(s,t);

      double temp = 1.
                  + 2. * alpha * zs * scale
                  + 2. * beta * scale * std::pow(zs, 1.5) * std::sin(3. * thetas)
                  + 2. * gamma * scale *  zs*zs
                  + 2. * delta * scale * std::pow(zs, 2.5) * std::sin(3.*thetas);
      return Norm * sqrt(xr * temp);
    };
};

#endif
