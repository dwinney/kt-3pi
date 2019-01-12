// Functions for Omnes function procedures to solve unitarity integral equations.
//
// Dependencies: pipi
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#include "omnes.hpp"

// Smoothly extrapolated Phase-shift shift matched at Lambda_phase^2
// Based off Danilkin et al, phi -> 3pi analysis
// http://cgl.soic.indiana.edu/jpac/w3pi.php
double omnes::extrap_phase(double s)
{
        double der, nJ, a, b;
        double LamSq = Lambda_phase*Lambda_phase;
        double match_pnt = phase_shift(wave, LamSq);

        double result;

        switch(wave) {
        case 0: {nJ = 2.; break;}
        case 1: {nJ = 1.; break;}
        };

        der = phase_shift(wave, LamSq + hD) - phase_shift(wave, LamSq - hD);
        der /= 2.*hD;

        a = 3. * pow(nJ*M_PI - match_pnt,2) / (2. * LamSq * der);
        b = 1. + 3. * (nJ*M_PI - match_pnt) / (2. * LamSq * der);

        if  (s <= LamSq)
        {
                result = phase_shift(wave, s);
        }
        else
        {
                result = nJ*M_PI - a / (b + pow(s/LamSq, 1.5));
        }

        return result;
};
