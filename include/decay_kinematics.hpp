// Generic kinematics class for the reaction of J^PC -> 3 pi.
// The reaction is defined by the mass and the quantum numbers of the decaying particle.
// This allows passing all relevant kinematic quantities easily into other classes.
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
#include <array>

using namespace std;

//-----------------------------------------------------------------------------
class decay_kinematics
{
protected:
  // Quantum Numbers
  int qn_J, qn_C, qn_P, qn_I, qn_H;

  double mDec;

  // Char string names of the decay particle (e.g. "omega")
  string decay_particle;

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
decay_particle(old.decay_particle)
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
};

//-----------------------------------------------------------------------------
// quantum_numbers.cpp

// Set functions
void set_decayJPC(int j, int p, int c);
void set_decayIsospin(int i);
void set_decayHelicity(int lambda);
void set_decayMass(double m);
void set_decayParticle(const char * n);

// Get functions
double get_decayMass();

int get_naturality();
int get_totalSpin();

array<int, 3> get_JPC();
string print_JPC();

string get_decayParticle();

//-----------------------------------------------------------------------------
// invariants.cpp

complex<double> t_man(complex<double> s, complex<double> zs); /// Mandelstam t in terms of s-channel CM variables
complex<double> u_man(complex<double> s, complex<double> t); // Mandelstam u

// Symmetric dimensionless variables for 3pi final state
double x(double s, double t);
double y(double s, double t);

// Lorentz Invariant Kibble Function
complex<double> Kibble(complex<double> s, complex<double> t);

// Dalitz region limits
double smin();
double smax();
double tmin(double s);
double tmax(double s);

// Geometric center of the Dalitz plot in s and t
double s_c();
double t_c();

//-----------------------------------------------------------------------------
// CoM_frame.cpp

// Triangle functions and aliases for
complex<double> Kallen(complex<double> x, complex<double> y, complex<double> z); // Kallen triangle function
complex<double> Kallen_x(complex<double> s);
complex<double> Kallen_pi(complex<double> s);

// Cosine of scattering angles in the s, t, and u center-of-mass frames
double z_s(double s, double t);
double z_t(double s, double t);
double z_u(double s, double t);

// Momenta in center of mass frame
complex<double> com_E2(double s);
complex<double> com_E3(double s);
complex<double> com_P2(double s);
complex<double> com_P3(double s);

// Center-of-mass scattering thresholds
double threshold();
double pseudo_threshold();

//-----------------------------------------------------------------------------
// angular_functions.cpp

// Outputs the Wigner little-d function appropriate for the decay into three scalar
double d_func0(int j, int l, double z);

// KSF Wigner-d functions
double d_hat(int j, int l, double z);
complex<double> d_hat(int j, int l, complex<double> z);

//-----------------------------------------------------------------------------
// kinematic_singularities.cpp

complex<double> barrier_factor(int j, int lam, complex<double> s);
complex<double> K_jlam(int j, int lam, complex<double> s, complex<double> zs);

//-----------------------------------------------------------------------------
// polar.cpp

double x_polar(double r, double phi);
double y_polar(double r, double phi);

double r(double s, double t);
double phi(double s, double t);
};

#endif
