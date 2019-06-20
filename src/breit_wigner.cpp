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
// KLOE BW amplitude
// ---------------------------------------------------------------------------
// Outputs the full amplitude, i.e. the sum of the pole in each sub-channel
// (NOTE: does not have the kinematic prefators associated with spin)
complex<double> breit_wigner::operator ()(double s, double t)
{
  double u = kinematics.u_man(s,t);
  complex<double> temp = F(s) + F(t) + F(u);
  return temp;
}

// Single-channel amplitude
complex<double> breit_wigner::F(double s)
{
  complex<double> denom = xr * (s - s_res()) + xi * sqrt(xr * s) * width(s);
  return s_res() / denom;
};

// Width with corrections from [https://doi.org/10.1103/PhysRevLett.21.244]
complex<double> breit_wigner::width(double s)
{
  complex<double> temp1, temp2;
  temp1 = mom_pi(s) / mom_pi(s_res());
  temp2 =   res_width * pow(temp1, 3.) * (s_res() / s);
  return temp2;
};

// Elastic pion momenta in the rho rest frame
complex<double> breit_wigner::mom_pi(double s)
{
   return .5 * sqrt(xr * (s - sthPi));
};

// Export a .dat file from threshold to 2 GeV
void breit_wigner::plot()
{
  std::ofstream output;

  std::string filename = "./BW_" + kinematics.get_ampName() + ".dat";
  output.open(filename.c_str());

  double s[100], re[100], im[100];

  double step = (2. - sthPi)/100.;
  complex<double> amp;

  for (int i = 0; i < 100; i++)
  {
    s[i] = sthPi + double(i) * step;
    amp = F(s[i]);
    re[i] = std::real(amp); im[i] = std::imag(amp);
    output << std::left << setw(15) << s[i] << setw(15) << re[i] << setw(15) << im[i] <<
    setw(15) << abs(amp * amp) << endl;
  }
output.close();
};

void breit_wigner::set_params(int n, const double * par)
{
  if (n > 2 || n < 0)
  {
    cout << "breit_wigner: Invalid number of free parameeters. Quitting... \n";
    exit(1);
  }
  switch(n)
  {
  case 2: {res_width = par[1];}
  case 1: {res_mass = par[0]; break;}
  }
};

void breit_wigner::print_params()
{
  if (kinematics.get_ampName() != "")
  {
    cout << kinematics.get_ampName() + ":";
  }
  cout << "Printing Breit-Wigner parameters... \n";
  cout << "---------------------------------------- \n";
  cout << std::left <<  setw(20) << "Mass:" << setw(20) << res_mass << "\n";
  cout << std::left <<  setw(20) << "Width:" << setw(20) << res_width << "\n";
  cout << "---------------------------------------- \n";
  cout << "\n";
};

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
  return f(s) + f(t) + f(kinematics.u_man(s,t));
}
