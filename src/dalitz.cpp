// General routines for Dalitz plot generation
//
// Dependencies: amp.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dalitz.hpp"

//-----------------------------------------------------------------------------
// Lorentz Invariant dimensionless parameters in terms of Mandelstam variables
template <class T>
double dalitz<T>::x(double s, double t)
{
  temp1 = sqrt(3.) * (t - u(s,t) );
  temp2 = 2.*mDec* (mDec - 3.*mPi);

  return temp1 / temp2;
};

template <class T>
double dalitz<T>::y(double s, double t)
{
  temp1 = 3. * (s_c() - s);
  temp2 = 2.* mDec * (mDec - 3.*mPi);

  return temp1 / temp2;
};

// Lorentz Invariant dimensionless parameters in terms of polar variables
template <class T>
double dalitz<T>::x_polar(double z, double theta)
{
  return sqrt(z) * cos(theta);
};

template <class T>
double dalitz<T>::y_polar(double z, double theta)
{
  return sqrt(z) * sin(theta);
};

// Inverted to get z and theta in terms of s and t
template <class T>
double dalitz<T>::z(double s, double t)
{
    temp1 = x(s,t);
    temp2 = y(s,t);

    return temp1 * temp1 + temp2 * temp2;
};

template <class T>
double dalitz<T>::theta(double s, double t)
{
  temp1 = x(s,t);
  temp2 = y(s,t);

  temp3 = atan2(temp2, temp1);
  if (temp3 < 0.0)
  {
    temp3 += 2. * M_PI; // Map [-pi,pi] to [0, 2pi]
  }
  return temp3;
};
//-----------------------------------------------------------------------------

// Doubly Differential Decay Width
template <class T>
double dalitz<T>::d2Gamma(double s, double t)
{
  temp1 = 96. * (mDec*mDec*mDec) * pow(2.* M_PI, 3);
  temp2 = abs(amp(s,t));
  return temp2 * temp2 / temp1;
};
