// Classes for General Amplitude Stuff including kinematic functions such as momenta.
//
// Dependencies: amp.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amp.hpp"

// Mandelstam u in terms of the other two variables
double amplitude::u(double s, double t)
{
 return 3.*mPi*mPi + mDec*mDec - s - t;
}

// Lorentz Covariant Kibble function, only defined in the physical regions, i.e. it can only be real
double amplitude::Kibble(double s, double t)
{
    double u = 3.*mPi*mPi + mDec*mDec - s - t;
    return (s * t * u ) - mPi*mPi * (mDec*mDec - mPi*mPi)*(mDec*mDec - mPi*mPi);
};

double amplitude::Kallen(double x, double y, double z)
{
  return x*x + y*y + z*z - 2.*(x*y + x*z + z*y);
};

// ---------------------------------------------------------------------------
// Momenta and Energies in the Center of Momentum Frame
double amplitude::com_E2(double s)
{
  return sqrt(s) / 2.;
};

double amplitude::com_E3(double s)
{
  return (mDec*mDec - mPi*mPi - s) / (2.*sqrt(s));
};

double amplitude::com_P2(double s)
{
  double temp = com_E2(s)*com_E2(s) - mPi*mPi;
  if (temp < 0.)
  {
    cout << " Outside of Dalitz region: com_P2 (sqrt(" << temp << ")) is imaginary. Quitting... \n";
    exit(1);
  }
  return sqrt(temp);
};

double amplitude::com_P3(double s)
{
  double temp = com_E3(s)*com_E3(s) - mPi*mPi;
  if (temp < 0.)
  {
    cout << " Outside of Dalitz region: com_P3 (sqrt(" << temp << ")) is imaginary. Quitting... \n";
    exit(1);
  }
  return sqrt(temp);
};

double amplitude::tmin(double s)
{
  return pow(com_E2(s) + com_E3(s), 2.) - pow(com_P2(s) + com_P3(s), 2.);
};

double amplitude::tmax(double s)
{
  return pow(com_E2(s) + com_E3(s), 2.) - pow(com_P2(s) - com_P3(s), 2.);
};
// ---------------------------------------------------------------------------
