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
public:
//-----------------------------------------------------------------------------
int qn_J, qn_C, qn_P, qn_I, qn_H;
double mDec;
double decay_particle;

void set_JPC(int j, int p, int c)
{
        qn_J = j; qn_P = p; qn_C = c;
};
void set_Isospin(int i)
{
        qn_I = i;
};
void set_Helicity(int lambda)
{
        qn_H = lambda;
}
void set_Mass(double m)
{
        mDec = m;
        cout << " Decay Mass set to: " << mDec << "\n";
        cout << endl;
};

//-----------------------------------------------------------------------------
// Kinematic Functions
//-----------------------------------------------------------------------------

double u(double s, double t); // Mandelstam u
double Kibble(double s, double t); // Lorentz Invariant Kibble Function
double Kallen(double x, double y, double z); // Kallen triangle function

double com_E2(double s);
double com_E3(double s);
double com_P2(double s);
double com_P3(double s);

// Dalitz region
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
//-----------------------------------------------------------------------------
};

#endif
