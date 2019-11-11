// Methods to evaluate invariant Mandelstam variables
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "decay_kinematics.hpp"

// ---------------------------------------------------------------------------
// Mandelstam t in terms of s-channel CM frame variables
complex<double> decay_kinematics::t_man(complex<double> s, complex<double> zs)
{
  complex<double> ks = sqrt(Kallen_pi(s) * Kallen_x(s)) / (4. * s);
  return (3. * mPi * mPi + mDec * mDec  - s) / 2. + 2. * ks * zs;
};

// Mandelstam u in terms of the other two invariant variables
complex<double> decay_kinematics::u_man(complex<double> s, complex<double> t)
{
 return 3.*mPi*mPi + mDec*mDec - s - t;
};

// ---------------------------------------------------------------------------
// Dalitz plot boundaries

// Lorentz Covariant Kibble function
complex<double> decay_kinematics::Kibble(complex<double> s, complex<double> t)
{
    complex<double> u = u_man(s,t);
    complex<double> lambda = (s * t * u) - mPi*mPi * (mDec*mDec - mPi*mPi)*(mDec*mDec - mPi*mPi);
    return xr * lambda;
};

double decay_kinematics::smin()
{
  return sthPi;
};

double decay_kinematics::smax()
{
  return (mDec - mPi) * (mDec - mPi);
};

double decay_kinematics::tmin(double s)
{
  return real(pow(com_E2(s) + com_E3(s), 2.) - pow(com_P2(s) + com_P3(s), 2.));
};

double decay_kinematics::tmax(double s)
{
  return real(pow(com_E2(s) + com_E3(s), 2.) - pow(com_P2(s) - com_P3(s), 2.));
};

// Center of the Dalitz plot in s and t
double decay_kinematics::s_c()
{
  return  (mDec * mDec + 3. * mPi* mPi) / 3.;
};

double decay_kinematics::t_c()
{
  return (3.*mPi*mPi + mDec*mDec - s_c())/2.;
};

// ---------------------------------------------------------------------------
// Lorentz Invariant dimensionless parameters in terms of Mandelstam variables

double decay_kinematics::x(double s, double t)
{
  double temp;
  temp = sqrt(3.) * real((t - u_man(s,t)));
  return temp / (2. * mDec * (mDec - 3. * mPi));
};

double decay_kinematics::y(double s, double t)
{
  double temp;
  temp = 3. * (s_c() - s);
  return temp / (2. * mDec * (mDec - 3. * mPi));
};
