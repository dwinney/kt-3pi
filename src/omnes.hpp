// Functions for Omnes function procedures to solve unitarity integral equations.
//
// Dependencies: pipi
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#ifndef _OMNES_
#define _OMNES_

#include "pipi.hpp"

class omnes : public pipi
{

protected:
int wave;

private:
static const double Lambda_phase = 1.3;
static const double hD = 0.0001;

public:
omnes() : pipi()
{

};
omnes(int i, int j) : pipi(i)
{
        wave = j;
};

double extrap_phase(double s);

};
#endif
