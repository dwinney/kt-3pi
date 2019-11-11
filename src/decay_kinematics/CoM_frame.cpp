// Methods to evaluate center-of-mass energies and momenta of each particle
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "decay_kinematics.hpp"

// ---------------------------------------------------------------------------
// Kallen functions
complex<double> decay_kinematics::Kallen(complex<double> x, complex<double> y, complex<double> z)
{
  return x*x + y*y + z*z - 2.*(x*y + x*z + z*y);
};

complex<double> decay_kinematics::Kallen_x(complex<double> s)
{
  return Kallen(s, mDec*mDec, mPi*mPi);
};

complex<double> decay_kinematics::Kallen_pi(complex<double> s)
{
  return Kallen(s, mPi*mPi, mPi*mPi);
};

//-----------------------------------------------------------------------------
// Scattering angles in the s, t, and u center-of-mass frames
double decay_kinematics::z_s(double s, double t)
{
  complex<double> ks = sqrt(Kallen_x(s) * Kallen_pi(s)) / s;
  complex<double> result =  (t - u_man(s,t)) / ks;

  if (abs(imag(result)) > 0.0001)
  {
    cout << "z_s: error arguments are complex. Quitting..." << endl;
    exit(1);
  }

  return real(result);
};

double decay_kinematics::z_t(double s, double t)
{
  return z_s(t , s);
};

double decay_kinematics::z_u(double s, double t)
{
  double u = real(u_man(s,t));
  return z_s(u, t);
};


// ---------------------------------------------------------------------------
// Momenta and Energies in the Center of Momentum Frame
complex<double> decay_kinematics::com_E2(double s)
{
  return sqrt(s * xr) / 2.;
};

complex<double> decay_kinematics::com_E3(double s)
{
  return (mDec*mDec - mPi*mPi - s * xr) / (2.*sqrt(s * xr));
};

complex<double> decay_kinematics::com_P2(double s)
{
  complex<double> temp = com_E2(s)*com_E2(s) - mPi*mPi;
  return sqrt(temp);
};

complex<double> decay_kinematics::com_P3(double s)
{
  complex<double> temp = com_E3(s)*com_E3(s) - mPi*mPi;
  return sqrt(temp);
};

// ---------------------------------------------------------------------------
// Regular initial state threshold
double decay_kinematics::threshold()
{
  return (mDec + mPi) * (mDec + mPi);
};

// Initial state pseudo-threshold
double decay_kinematics::pseudo_threshold()
{
  return  (mDec - mPi) * (mDec - mPi);
};
