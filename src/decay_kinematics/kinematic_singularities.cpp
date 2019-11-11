// Methods to evaluate the kinematic singularitity contributions of a given amplitude
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "decay_kinematics.hpp"

//-----------------------------------------------------------------------------
// Angular momentum barrier factor
complex<double> decay_kinematics::barrier_factor(int j, int lam, complex<double> s)
{
  if (j < std::abs(lam))
  {
    cout << "K_jlam: j < |lambda| is not allowed. Quitting..." << endl;
    exit(1);
  }

  double dlam = double(lam), dj = double(j);

  // Angular momentum barrier factor
  complex<double> temp2 = sqrt(Kallen_x(s) * Kallen_pi(s));
  return pow(temp2, dj - abs(dlam));
};

//-----------------------------------------------------------------------------
// Kinematic Singularities of Helicity amplitudes
// helicity projection : Lambda
// spin projection : j > 0
complex<double> decay_kinematics::K_jlam(int j, int lam, complex<double> s, complex<double> zs)
{
  if (j < std::abs(lam))
  {
    cout << "K_jlam: j < |lambda| is not allowed. Quitting..." << endl;
    exit(1);
  }

  double dlam = double(lam), dj = double(j);

  // Half angle factors and K factors
  complex<double> temp = Kibble(s, t_man(s, zs)) / 4.;
  complex<double> result = pow(sqrt(temp), abs(dlam));

  // LS - lambda little group factor
  int Yx = qn_J - (1 + get_naturality())/2;
  complex<double> temp3 = sqrt(Kallen_x(s));
  result *= pow(temp3, abs(dj - double(Yx)) - dj);

  return result;
};
