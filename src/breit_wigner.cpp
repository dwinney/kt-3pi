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
complex<double> breit_wigner_KLOE::operator ()(double s, double t)
{
  double real_part, imag_part;

  real_part = s - s_res();
  imag_part = sqrt(s) * width(s);

  return s_res() / (real_part * xr + imag_part * xi);
};

double breit_wigner_KLOE::width(double s)
{
  double temp1, temp2;

  temp1 = com_P3(s) / com_P3(s_res());
  temp2 = res_width * pow(temp1, 3.) * (s_res() / s);

  return temp2;
};
