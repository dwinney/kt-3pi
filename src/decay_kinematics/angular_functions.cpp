// Rotational functions for evaluating amplitudes in different frames
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "decay_kinematics.hpp"

//-----------------------------------------------------------------------------
// Angular functions

// Outputs the Wigner little-d function appropriate for the decay into three scalar,
// d^j_{lambda, 0}(z)
double decay_kinematics::d_func0(int j, int l, double z)
{
  double prefactor = TMath::Factorial(j - l) / TMath::Factorial(j + l);

  double legendre;
  if (l >= 0)
  {
    legendre = ROOT::Math::assoc_legendre(j, l, z);
  }
  else
  {
    legendre = ROOT::Math::assoc_legendre(j, -l, z);
    legendre *= pow(-1., double(l));
    legendre *= TMath::Factorial(j - l) / TMath::Factorial(j + l);
  }

  return sqrt(2.) * sqrt(prefactor) * legendre;
};

// Outputs d_hat the kinematic-singularity-free d_function
double decay_kinematics::d_hat(int j, int l, double z)
{
  double d_lam = d_func0(j, l, z);
  double half_angle_factor = pow((1. - z*z), abs(double(l))/2.);

  return d_lam / half_angle_factor;
};

complex<double> decay_kinematics::d_hat(int j, int l, complex<double> z)
{
  if (l == 1)
  {
    switch (j)
    {
    case 1: return 1.;
    case 3: return 3. * z;
    default: cout << "d_hat: complex legendres not coded that far yet... \n";
              exit(0);
    }
  }
  else
  {
    cout << "d_hat: complex legendres not coded that far yet... \n";
    exit(0);
  }
};
