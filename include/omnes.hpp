// Functions for Omnes function procedures.
//
// Dependencies: pipi, aux_math, amp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Define an omnes object with a given partial wave and isospin:
// omnes my_omnes(isospin, spin);
//
// Evaluate the Omnes function (as a complex number) at some s:
// my_omnes(s);
// -----------------------------------------------------------------------------

#ifndef _OMNES_
#define _OMNES_

#include "pipi.hpp"
#include "decay_kinematics.hpp"
#include "aux_math.hpp"

//-----------------------------------------------------------------------------
class omnes : public pipi, public decay_kinematics
//-----------------------------------------------------------------------------
{
protected:
int wave;
int sign = 1;
double eps = 1e-9;
int N_omnes = 60;
double s0 = sthPi; //Lower Bound for integral

//-----------------------------------------------------------------------------
private:
static constexpr double Lambda_phase = 1.3;
static constexpr double LamSq = Lambda_phase*Lambda_phase;
static constexpr double hD = 0.0001;

bool WG_GENERATED;
std::vector<double> wgt, abs;

//-----------------------------------------------------------------------------
public:
// Default constructor
omnes() : pipi()
{
};

// Parameterized constructor with quantum numbers
omnes(int i, int j, const char * n = "") : pipi(i), wave(j)
{
  set_ampName(n);
};

// Copy constructor
omnes(const omnes &previous):
  pipi(previous.pipi::pipi_qn_I), wave(previous.wave),
  wgt(previous.wgt), abs(previous.abs), WG_GENERATED(previous.WG_GENERATED)
{};

std::complex<double> operator ()(double s, double t);

void set_ieps(int sig);
void set_N_omnes(int i);

double extrap_phase(double s);
double kernel(double s, double sp);
std::complex<double> eval(double s);

};
//-----------------------------------------------------------------------------
#endif
