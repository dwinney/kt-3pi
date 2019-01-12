// Self-contained implementation of the GKPY paramterization for low-energy pion scattering.
// Based on: 10.1103/PhysRevD.83.074004
//
// Dependencies: None
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _PIPI_
#define _PIPI_

#include <complex>
#include <iostream>
#include <stdlib.h>
#include <cmath>

//-----------------------------------------------------------------------------
// Generic Amplitude Object to store relevant Quantum Numbers
class amplitude
{
protected:
int qn_J, qn_C, qn_P, qn_I, qn_H;     // quantum numbers

public:
void setJPC(int j, int p, int c)
{
        qn_J = j; qn_P = p; qn_C = c;
};
void setIsospin(int i)
{
        qn_I = i;
};
void setHelicity(int lambda)
{
        qn_H = lambda;
}
};
//-----------------------------------------------------------------------------

// Pion Scattering Amplitude Object
class pipi : public amplitude
{

public:
static const double conv; // degrees to radians conversion
static const std::complex<double> xr, xi; // unit imaginary and real

static const double mPi, mK, mEta, mRho, mF2; //masses
static const double sthPi, sthK, sthEta; // particle thresholds

pipi();
pipi(int i);

double conformal(double s, double s0);
double elastic_mom( double s, double sth);
double legendre(int l, double x);

double phase_shift(int l, double s);
double inelasticity(int l, double s);

std::complex<double> GKPRY_partial_wave(int l, double s);
std::complex<double> GKPRY_iso_amp(double s, double z);
};

#endif
