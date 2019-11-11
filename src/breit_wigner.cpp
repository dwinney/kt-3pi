// Class for different Breit-Wigner parameterizations of amplitudes
//
// Dependencies: decay_kinematics.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "breit_wigner.hpp"

// Outputs the full amplitude, i.e. the sum of the pole in each sub-channel
// (NOTE: does not have the kinematic prefators associated with spin)
complex<double> breit_wigner::eval(double s, double t)
{

  double u = real(kinematics.u_man(s,t));
  complex<double> temp = F(s) + F(t) + F(u);

  return kinematics.K_jlam(1, 1, s, kinematics.z_s(s,t)) * temp;
}

// Single-channel amplitude
complex<double> breit_wigner::F(double s)
{
  complex<double> denom = xr * (s - s_res()) + xi * sqrt(xr * s) * width(s);
  return normalization * s_res() / denom;
};

// Width with corrections from
// [https://doi.org/10.1103/PhysRevLett.21.244]
complex<double> breit_wigner::width(double s)
{
  complex<double> temp1, temp2;
  temp1 = mom_pi(s) / mom_pi(s_res());
  temp2 = res_width * pow(temp1, 3.) * (s_res() / s);
  return temp2;
};

// Elastic pion momenta in the rho rest frame
complex<double> breit_wigner::mom_pi(double s)
{
   return .5 * sqrt(xr * (s - sthPi));
};

// Export a .dat file from threshold to 1 GeV
void breit_wigner::plot()
{
  std::ofstream output;

  std::string filename = "./BW.dat";
  output.open(filename.c_str());

  std::vector<double> s;
  std::vector<complex<double>> fx;

  double step = (1. - sthPi)/100.;

  for (int i = 0; i < 100; i++)
  {
    double s_i = sthPi + double(i) * step;
    s.push_back(sqrt(s_i));

    complex<double> amp = - F(s_i);
    fx.push_back(amp);

    output << std::left << setw(15) << sqrt(s_i) << setw(15) << std::real(amp) << setw(15) << std::imag(amp) << endl;
  }

  output.close();

  quick_plot(s, fx, "BW");
};

// ----------------------------------------------------------------------------
// Set the normalization coefficient
void breit_wigner::normalize(double gamma_exp)
{
  cout << endl;
  cout << "Normalizing Breit-Wigner amplitude to Gamma_3pi = " << gamma_exp << " MeV..." << endl;

  dalitz d_plot(this);
  double gamma = d_plot.decay_width();
  normalization = sqrt(gamma_exp * 1.e-3 / gamma);

  cout << "Normalization constant = " << normalization << "." << endl;
  cout << endl;
};

void breit_wigner::set_params(int n, const double * par)
{
  if (n > 2 || n < 0)
  {
    cout << "breit_wigner: Invalid number of free parameters. Quitting... \n";
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
  cout << "Printing Breit-Wigner parameters... \n";
  cout << "---------------------------------------- \n";
  cout << std::left <<  setw(20) << "Mass:" << setw(20) << res_mass << "\n";
  cout << std::left <<  setw(20) << "Width:" << setw(20) << res_width << "\n";
  cout << "---------------------------------------- \n";
  cout << "\n";
};
