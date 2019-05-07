// Classes for Dalitz plot generation and fitting to polynomial expansion with dalitz plot parameters.
//
//
// Dependencies: amp.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dalitz.hpp"

//-----------------------------------------------------------------------------
// Lorentz Invariant dimensionless parameters in terms of Mandelstam variables
double dalitz::x(double s, double t)
{
  temp_num = std::sqrt(3.) * (t - u(s,t) );
  temp_den = 2.*Mass* (Mass - 3.*mPi);

  return temp_num / temp_den;
};

double dalitz::y(double s, double t)
{
  temp_num = 3. * (s_c - s);
  temp_den = 2.* Mass * (Mass - 3.*mPi);
};

// Lorentz Invariant dimensionless parameters in terms of polar variables
double dalitz::x_polar(double z, double theta)
{
  return std::sqrt(z) * std::cos(theta);
};

double dalitz::y_polar(double z, double theta)
{
  return std::sqrt(z) * std::sin(theta);
};
//-----------------------------------------------------------------------------

// Lorentz Invariant dimensionless parameters in terms of polar variables
double dalitz::d2Gamma(double s, double t)
{
  temp_pref = 96. * Mass*Mass*Mass * std::pow(2.*M_PI, 3);
  amp_sqr = real(amp(s,t) * conj(amp(s,t)));
  return temp_pref;
};
