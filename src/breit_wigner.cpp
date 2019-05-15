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
complex<double> breit_wigner_simple::f(double s)
{
  double real_part, imag_part;
  real_part = s - s_res();
  imag_part = res_width * res_mass;

  return 1. / (real_part * xr - imag_part * xi);
};

complex<double> breit_wigner_simple::operator ()(double s, double t)
{
  return f(s) + f(t) + f(u(s,t));
}

// ---------------------------------------------------------------------------
// KLOE BW amplitude
// ---------------------------------------------------------------------------
complex<double> breit_wigner_KLOE::f(double s)
{
  complex<double> denom = s - s_res() + xi * res_mass * sqrt(s) * width(s);
  return s_res() / denom;
};

complex<double> breit_wigner_KLOE::operator ()(double s, double t)
{
  double u_man = amplitude::u(s,t);
  complex<double> temp = f(s) + f(t) + f(u_man);
  return temp;
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
   return .5 * sqrt( xr * (s - sthPi));
};

void breit_wigner_KLOE::plot()
{
  std::ofstream output;
  std::string filename = "./BW_KLOE_" + name + ".dat";
  output.open(filename.c_str());
  double s[100], re[100], im[100];

  double step = (2. - sthPi)/100.;
  complex<double> amp;
  for (int i = 0; i < 100; i++)
  {
    s[i] = sthPi + double(i) * step;
    amp = f(s[i]);
    re[i] = std::real(amp); im[i] = std::imag(amp);
    output << std::left << setw(15) << s[i] << setw(15) << re[i] << setw(15) << im[i] <<
    setw(15) << abs(amp * amp) << endl;
  }
output.close();
};
