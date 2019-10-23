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

#include <TMath.h>
#include <Math/SpecFuncMathMore.h>

#include <iostream>
#include <complex>
#include <string>

using std::cout;
using std::endl;
using std::string;

//-----------------------------------------------------------------------------
class decay_kinematics
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

int get_naturality()
{
  if (qn_P == pow(-1, qn_J))
  {
    return 1;
  }
  else
  {
    return -1;
  }
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

  if (qn_P == 1)
  {
    result += "+";
  }
  else if (qn_P == -1)
  {
    result += "-";
  }

  if (qn_C == 1)
  {
    result += "+";
  }
  else if (qn_C == -1)
  {
    result += "-";
  }

  return result;
};

int get_totalSpin()
{
  return qn_J;
}

//-----------------------------------------------------------------------------
// Kinematic Functions
complex<double> t_man(complex<double> s, complex<double> zs); /// Mandelstam t in terms of s-channel CM variables
complex<double> u_man(complex<double> s, complex<double> t); // Mandelstam u

complex<double> Kibble(complex<double> s, complex<double> t); // Lorentz Invariant Kibble Function
complex<double> Kallen(complex<double> x, complex<double> y, complex<double> z); // Kallen triangle function

// aliases for the two relevant kallen functions
complex<double> Kallen_x(complex<double> s)
{
  return Kallen(s, mDec*mDec, mPi*mPi);
};

complex<double> Kallen_pi(complex<double> s)
{
  return Kallen(s, mPi*mPi, mPi*mPi);
};

//-----------------------------------------------------------------------------
// Cosine of scattering angles in the s, t, and u center-of-mass frames
double z_s(double s, double t);
double z_t(double s, double t);
double z_u(double s, double t);

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
// Angular functions

// Outputs the Wigner little-d function appropriate for the decay into three scalar,
// d^j_{lambda, 0}(z)
double d_func0(int j, int l, double z);

// Outputs d_hat the kinematic-singularity-free d_function
double d_hat(int j, int l, double z);
complex<double> d_hat(int j, int l, complex<double> z);

//-----------------------------------------------------------------------------
// Kinematic Singularities of Helicity amplitudes
complex<double> barrier_factor(int j, int lam, complex<double> s);
complex<double> K_jlam(int j, int lam, complex<double> s, complex<double> zs);

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
    return (3.*mPi*mPi + mDec*mDec - s_c())/2.;
};

//-----------------------------------------------------------------------------
// Lorentz Invariant dimensionless variables X and Y
double d_norm(){
  return 1. / (2. * mDec * (mDec - 3. * mPi));
}
double x(double s, double t);
double y(double s, double t);

double x_polar(double r, double phi);
double y_polar(double r, double phi);

double r(double s, double t);
double phi(double s, double t);
};

#endif
