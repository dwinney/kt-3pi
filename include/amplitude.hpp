// Generic amplitude class for the reaction of J^PC -> 3 pi.
// The reaction is defined by the mass and the quantum numbers of the decaying particle.
//
// Dependencies: constants.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _AMP_
#define _AMP_

#include "constants.hpp"
#include <iostream>
#include <complex>
#include <string>

using std::complex;
using std::cout;
using std::endl;
using std::string;

//-----------------------------------------------------------------------------
class amplitude
//-----------------------------------------------------------------------------
{
protected:
  int qn_J, qn_C, qn_P, qn_I, qn_H;
  double mDec;
  string decay_particle, amp_name;

public:
//-----------------------------------------------------------------------------
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
        // 
        // if (amp_name != "")
        // {
        //   cout << amp_name + ": ";
        // }
        // cout << "Decay Mass set to: " << mDec << "\n";
        // cout << endl;
};

void set_ampName(const char * n)
{
  amp_name = n;
};

void set_decayParticle(const char * n)
{
  decay_particle = n;
};

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

//-----------------------------------------------------------------------------
// Kinematic Functions
//-----------------------------------------------------------------------------

double u_man(double s, double t); // Mandelstam u
double Kallen(double x, double y, double z); // Kallen triangle function

complex<double> Kibble(double s, double t); // Lorentz Invariant Kibble Function

// Momenta in center of mass frame
complex<double> com_E2(double s);
complex<double> com_E3(double s);
complex<double> com_P2(double s);
complex<double> com_P3(double s);

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
// Lorentz Invariant dimensionless parameters X and Y
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
