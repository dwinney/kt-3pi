// Functions for Omnes function procedures.
//
// Dependencies: pipi, aux_math
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#ifndef _OMNES_
#define _OMNES_

#include "pipi.hpp"
#include "aux_math.hpp"

//-----------------------------------------------------------------------------
class omnes : public pipi
//-----------------------------------------------------------------------------
{
protected:
int wave;
int sign = 1;
double eps = 1e-9;
int N_omnes = 100;
double s0 = sthPi; //Lower Bound for integral

//-----------------------------------------------------------------------------
private:
static constexpr double Lambda_phase = 1.3;
static constexpr double LamSq = Lambda_phase*Lambda_phase;
static constexpr double hD = 0.0001;

bool WG_GENERATED = false;
std::vector<double> wgt, abs;

//-----------------------------------------------------------------------------
public:
omnes() : pipi()
{
};
omnes(int i, int j) : pipi(i)
{
        wave = j;
};

void set_eps(int sig, double eps);
void set_N_omnes(int i);

double extrap_phase(double s);
double kernel(double s, double sp);
std::complex<double> eval(double s);

};
//-----------------------------------------------------------------------------
#endif
