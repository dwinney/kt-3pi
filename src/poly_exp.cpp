// Simple class with a polynomial amplitude mimicing the expansion around center of dalitz plot to test fitting and parameter extraction.
//
// Dependencies: decay_kinematics
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "poly_exp.hpp"

complex<double> poly_exp::eval(double s, double t)
{
  double rs = kinematics.r(s,t);
  double phis =  kinematics.phi(s,t);

  // cout << setw(10) << s << setw(10) << t << setw(10) << rs << endl;
  complex<double> temp = xr;
  temp += 2. * alpha * rs * scale;
  temp += 2. * beta * scale * pow(rs, 1.5) * sin(3. * phis);
  temp += 2. * gamma * scale * rs * rs;
  temp += 2. * delta * scale * pow(rs, 2.5) * std::sin(3. * phis);

  complex<double> result = Norm * sqrt(temp);
  result *= kinematics.K_jlam(1, 1, s, kinematics.z_s(s,t));

  return result;
};

void poly_exp::set_params(int n, const double *par)
{
  switch (n - 1)
   {
    case 4: delta = par[4];
    case 3: gamma = par[3];
    case 2: beta = par[2];
    case 1: alpha = par[1];
    case 0: Norm = par[0]; break;
  };
};

void poly_exp::set_errors(int n, const double *par)
{
  switch (n - 1)
   {
    case 4: ddelta = par[4];
    case 3: dgamma = par[3];
    case 2: dbeta = par[2];
    case 1: dalpha = par[1];
    case 0: dNorm = par[0]; break;
  };
};

double poly_exp::error_func(double s, double t)
{
  if (ERROR == true)
  {
    double result;
    double rs = kinematics.r(s,t);
    double phis =  kinematics.phi(s,t);

    result = pow(2. * rs * dalpha, 2.);
    result += rs * rs * rs * pow(2. * std::sin(3. * phis) * dbeta, 2.);
    result += pow(2. * rs * rs * dgamma, 2.);

    return sqrt(result);
  }
  else
  {
    return 1.;
  };
};

void poly_exp::print_params(int a)
{
  switch (a)
  {
  case 0: {
      cout << "Printing Dalitz Plot parameters... \n";
      cout << "---------------------------------------- \n";
      cout << std::left <<  setw(20) << "normalization:" << setw(20) << Norm << "\n";
      cout << std::left <<  setw(20) << "alpha:" << setw(20) << alpha << "\n";
      cout << std::left <<  setw(20) << "beta:" << setw(20) << beta << "\n";
      cout << std::left <<  setw(20) << "gamma:" << setw(20) <<  gamma << "\n";
      cout << std::left <<  setw(20) << "delta:" << setw(20) << delta << "\n";
      cout << "---------------------------------------- \n";
      cout << "\n";
    };
  };
};

void poly_exp::normalize(double gamma_exp)
{
  dalitz d_plot(this);
  double gamma = d_plot.Gamma_total();
  Norm = sqrt(gamma_exp * 1.e-3 / gamma);
};
