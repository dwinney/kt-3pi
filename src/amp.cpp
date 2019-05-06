// Classes for General Amplitude Stuff
//
// Dependencies: amp.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amp.hpp"

//-----------------------------------------------------------------------------
double amplitude::conv = (M_PI / 180.);

//Masses
double amplitude::mPi = 0.1396;
double amplitude::mK = 0.496;
double amplitude::mEta = 0.54753;

double amplitude::mRho = .77545;
double amplitude::mF2 = 1.2754;

//Thresholds for pi, eta, and K
double amplitude::sthPi = 4.*mPi*mPi;
double amplitude::sthK = 4.*mK*mK;
double amplitude::sthEta = 4.*mEta*mEta;

//Unit imaginary and real
std::complex<double> amplitude::xr(1., 0.);
std::complex<double> amplitude::xi(0., 1.);

// Lorentz Covariant Kibble function, only defined in the physical reqgions, i.e. it can only be real
double amplitude::Kibble(double s, double t)
{
 return 1;
};
