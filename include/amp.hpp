// Classes for General Amplitude Stuff
//
// Dependencies: None
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

//-----------------------------------------------------------------------------
// Generic Amplitude Object to store relevant Quantum Numbers for the Decaying Particle
class amplitude
//-----------------------------------------------------------------------------
{
protected:
int qn_J, qn_C, qn_P, qn_I, qn_H;
double Mass;
//-----------------------------------------------------------------------------
public:
amplitude(){};

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
        Mass = m;
};

// Kinematic Functions
double u(double s, double t);
double Kibble(double s, double t);
//-----------------------------------------------------------------------------
};

#endif
