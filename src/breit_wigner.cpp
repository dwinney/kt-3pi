// Class for different Breit-Wigner parameterizations of amplitudes
//
// Dependencies: amp.cpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "breit_wigner.hpp"
// ---------------------------------------------------------------------------
// Simple Relativisitic Breit-Wigner amplitude with constant width
complex<double> breit_wigner_simple::operator ()(double s, double t)
{
  double real_part, imag_part;
  real_part = s - s_res();
  imag_part = res_width * res_mass;

  return 1. / (real_part * xr - imag_part * xi);
};

// ---------------------------------------------------------------------------
// KLOE BW amplitude
// ---------------------------------------------------------------------------
complex<double> breit_wigner_KLOE::f(double s)
{
  complex<double> denom = s - s_res() + xi * std::sqrt(xr * s) * width(s);
  return s_res() / denom;
};

complex<double> breit_wigner_KLOE::operator ()(double s, double t)
{
  double u_man = amplitude::u(s,t);
  return f(s) + f(t) + f(u_man);
}

complex<double> breit_wigner_KLOE::width(double s)
{
  complex<double> temp1, temp2;
  temp1 = mom_pi(s) / mom_pi(s_res());
  temp2 = res_width * pow(temp1, 3.) * (s_res() / s);
  return temp2;
};

complex<double> breit_wigner_KLOE::mom_pi(double s)
{
  complex<double> temp1, temp2;
   temp1 = (mDec*mDec - mPi*mPi - s) / (2.*std::sqrt(s * xr));
   temp2 = sqrt(temp1 * temp1 - mPi * mPi);

   return temp2;
};
