// Header file with global phyiscal constants. Everything is in GeV unless explicitly stated otherwise.
//
// Dependencies: None
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _CONSTANT_
#define _CONSTANT_

#include <cmath>
#include <complex>

using std::complex;

//-----------------------------------------------------------------------------
const double conv = (M_PI / 180.);
const double elastic_cutoff = 0.95;

//Masses
const double mPi = 0.1396;
const double mK = 0.496;
const double mEta = 0.54753;

const double mRho = .77545;
const double mF2 = 1.2754;

//Thresholds for pi, eta, and K
const double sthPi = 4.*mPi*mPi;
const double sthK = 4.*mK*mK;
const double sthEta = 4.*mEta*mEta;

//Unit imaginary and real
const complex<double> xr(1., 0.);
const complex<double> xi(0., 1.);

#endif
