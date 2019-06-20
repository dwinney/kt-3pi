// Simple class with a polynomial amplitude mimicing the expansion around center of dalitz plot to test fitting and parameter extraction.
//
// Dependencies: decay_kinematics
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "poly_exp.hpp"

complex<double> poly_exp::operator ()(double s, double t)
  {
    double zs = kinematics.z(s,t);
    double thetas =  kinematics.theta(s,t);

    complex<double> temp = 1.
                + 2. * alpha * zs * scale
                + 2. * beta * scale * std::pow(zs, 1.5) * std::sin(3. * thetas)
                + 2. * gamma * scale *  zs*zs
                + 2. * delta * scale * std::pow(zs, 2.5) * std::sin(3.*thetas);
    return Norm * sqrt(temp);
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
