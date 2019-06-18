// Generic kinematics class for the reaction of J^PC -> 3 pi.
// The reaction is defined by the mass and the quantum numbers of the decaying particle.
// This allows passing all relevant kinematic quantities easily into other classes.
//
// Dependencies: constants.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _DECAY_
#define _DECAY_

#include "constants.hpp"

#include <iostream>
#include <complex>
#include <string>

using std::cout;
using std::endl;
using std::string;

//-----------------------------------------------------------------------------
class decay_kinematics
//-----------------------------------------------------------------------------
{
protected:
  // Quantum Numbers
  int qn_J, qn_C, qn_P, qn_I, qn_H;
  double mDec;

  // Char string names of the decay particle (e.g. "omega") and
  // for the user-defined amplitude (e.g. "eff_breit_wigner")
  string decay_particle, amp_name;

public:
//-----------------------------------------------------------------------------
// Default Constructor
// For clarity and keeping track of which numbers go where in a script
// instead of having a parameterized constructor,
// I prefer the data members be set with the set and get functions.
decay_kinematics()
{};

// Copy Constructor
decay_kinematics(const decay_kinematics & old) :
qn_J(old.qn_J), qn_C(old.qn_C), qn_P(old.qn_P),
qn_H(old.qn_H), qn_I(old.qn_I),
mDec(old.mDec),
 decay_particle(old.decay_particle), amp_name(old.amp_name)
{};

// Equality operator to make copies as well
void operator = (const decay_kinematics &old)
{
  qn_J = old.qn_J;
  qn_C = old.qn_C;
  qn_P = old.qn_P;
  qn_H = old.qn_H;
  qn_I = old.qn_I;
  mDec = old.mDec;
  decay_particle = old.decay_particle;
  amp_name = old.amp_name;
}
void set_decayJPC(int j, int p, int c)
{
        qn_J = j; qn_P = p; qn_C = c;
};
void set_decayIsospin(int i)
{
        qn_I = i;
};
void set_decayHelicity(int lambda)
{
        qn_H = lambda;
}
void set_decayMass(double m)
{
        mDec = m;
};

void set_ampName(const char * n)
{
  amp_name = n;
};

void set_decayParticle(const char * n)
{
  decay_particle = n;
};

// Generic implementation-dependent function to set n free parameters
void set_params(int n, const double *par);

//-----------------------------------------------------------------------------
// get functions to access protected data outside of the initial construction
double get_decayMass()
{
  return mDec;
};

string get_ampName()
{
  return amp_name;
};

string get_decayParticle()
{
  return decay_particle;
};

string get_JPC()
{
  string result = std::to_string(qn_J);

  if (qn_P = 1)
  {
    result += "+";
  }
  else if (qn_P = -1)
  {
    result += "-";
  }

  if (qn_C = 1)
  {
    result += "+";
  }
  else if (qn_C = -1)
  {
    result += "-";
  }

  return result;
};
//-----------------------------------------------------------------------------
// Kinematic Functions
double u_man(double s, double t); // Mandelstam u
double Kallen(double x, double y, double z); // Kallen triangle function

complex<double> Kibble(double s, double t); // Lorentz Invariant Kibble Function

// Momenta in center of mass frame
complex<double> com_E2(double s);
complex<double> com_E3(double s);
complex<double> com_P2(double s);
complex<double> com_P3(double s);

// Threshold
double threshold(){
  return (mDec + mPi) * (mDec + mPi);
};

double pseudo_threshold()
{
  return  (mDec - mPi) * (mDec - mPi);
};

//-----------------------------------------------------------------------------
// Dalitz region limits
double smin()
  {
    return sthPi;
  };
double smax()
  {
    return (mDec - mPi) * (mDec - mPi);
  };
double tmin(double s);
double tmax(double s);

// Center of the Dalitz plot in s and t
double s_c()
  {
    return  (mDec * mDec + 3. * mPi* mPi) / 3.;
  };
double t_c()
{
    return (3.*mPi*mPi + mDec* mDec - s_c())/2.;
};

//-----------------------------------------------------------------------------
// Lorentz Invariant dimensionless variables X and Y
double d_norm(){
  return 1. / (2. * mDec * (mDec - 3. * mPi));
}
double x(double s, double t);
double y(double s, double t);

double x_polar(double z, double theta);
double y_polar(double z, double theta);

double z(double s, double t);
double theta(double s, double t);
};

#endif
