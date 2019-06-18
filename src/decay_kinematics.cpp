// Classes for General Amplitude Stuff including kinematic functions such as momenta.
//
// Dependencies: amp.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "decay_kinematics.hpp"

// Mandelstam u in terms of the other two variables
double decay_kinematics::u_man(double s, double t)
{
 return 3.*mPi*mPi + mDec*mDec - s - t;
};
// Lorentz Covariant Kibble function
complex<double> decay_kinematics::Kibble(double s, double t)
{
    double u = u_man(s,t);
    double lambda = (s * t * u ) - mPi*mPi * (mDec*mDec - mPi*mPi)*(mDec*mDec - mPi*mPi);
    return xr * lambda;
};

double decay_kinematics::Kallen(double x, double y, double z)
{
  return x*x + y*y + z*z - 2.*(x*y + x*z + z*y);
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

double decay_kinematics::tmin(double s)
{
  return real(pow(com_E2(s) + com_E3(s), 2.) - pow(com_P2(s) + com_P3(s), 2.));
};

double decay_kinematics::tmax(double s)
{
  return real(pow(com_E2(s) + com_E3(s), 2.) - pow(com_P2(s) - com_P3(s), 2.));
};
// ---------------------------------------------------------------------------
// Lorentz Invariant dimensionless parameters in terms of Mandelstam variables
double decay_kinematics::x(double s, double t)
{
  double temp;
  temp = sqrt(3.) * (t - u_man(s,t));
  return temp * d_norm();
};

double decay_kinematics::y(double s, double t)
{
  double temp;
  temp = 3. * (s_c() - s);
  return temp * d_norm();
};
//-----------------------------------------------------------------------------
// Lorentz Invariant dimensionless parameters in terms of polar variables
double decay_kinematics::x_polar(double z, double theta)
{
  return sqrt(z) * cos(theta);
};

double decay_kinematics::y_polar(double z, double theta)
{
  return sqrt(z) * sin(theta);
};
//-----------------------------------------------------------------------------
// Inverted to get z and theta in terms of s and t
double decay_kinematics::z(double s, double t)
{
    double tmp1, tmp2;
    tmp1 = x(s,t);
    tmp2 = y(s,t);
    return tmp1 * tmp1 + tmp2 * tmp2;
};

double decay_kinematics::theta(double s, double t)
{
  double temp1, temp2, temp3;
  temp1 = x(s,t);
  temp2 = y(s,t);

  temp3 = atan2(temp2, temp1);
  if (temp3 < 0.0)
  {
    temp3 += 2. * M_PI; // Map [-pi,pi] to [0, 2pi]
  }
  return temp3;
};