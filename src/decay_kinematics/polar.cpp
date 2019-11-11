// Invariant Mandelstam variables in terms of polar variables
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "decay_kinematics.hpp"

//-----------------------------------------------------------------------------
// Lorentz Invariant dimensionless parameters in terms of polar variables
double decay_kinematics::x_polar(double r, double phi)
{
  return sqrt(r) * cos(phi);
};

double decay_kinematics::y_polar(double r, double phi)
{
  return sqrt(r) * sin(phi);
};

//-----------------------------------------------------------------------------
// Inverted to get r and phi in terms of s and t
double decay_kinematics::r(double s, double t)
{
    double X, Y;
    X = x(s,t);
    Y = y(s,t);
    return X * X + Y * Y;
};

double decay_kinematics::phi(double s, double t)
{
  double result = atan2(y(s,t), x(s,t));
  if (result < 0.0)
  {
    result += 2. * M_PI; // Map [-pi,pi] to [0, 2pi]
  }
  return result;
};
