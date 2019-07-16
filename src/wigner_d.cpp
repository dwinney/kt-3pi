// Functions related to the Wigner little-d functions and the half-angle factors
// which remove the kinematic singularities in the cross-variable
//
// Dependencies: None
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#include "wigner_d.hpp"

double d_func(int j, int l, double z)
{
  double prefactor = TMath::Factorial(j - l) / TMath::Factorial(j + l);
  double legendre = ROOT::Math::assoc_legendre(j, l, z);

  return sqrt(prefactor) * legendre;
};

double d_hat(int j, int l, double z)
{
  double d_lam = d_func(j, l, z);
  double half_angle_factor = pow((1. - z*z), double(l)/2.);

  return d_lam / half_angle_factor;
};
