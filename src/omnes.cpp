// Functions for Omnes function procedures to solve unitarity integral equations.
//
// Dependencies: pipi
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#include "omnes.hpp"

// Set the i epsilon persription for the integration kernal in Omnes function
// Default is + 1e-9
void omnes::set_eps(int sig, double e)
{
        sign = sig;
        if (sig != 1 && sig !=-1 && sig != 0)
        {
                std::cout << " Invalid i*eps presription. Allowed sig values = -1, 0, +1"
                          << "\n";
                exit(1);
        }
        eps = e;
};

// Set number of Gaussian Quandrature points in evaluating Omnes function
void omnes::set_N_omnes(int i)
{
        N_omnes = i;
};

// Smoothly extrapolated Phase-shift shift matched at Lambda_phase^2
// Based off Danilkin et al, phi -> 3pi analysis
// http://cgl.soic.indiana.edu/jpac/w3pi.php
double omnes::extrap_phase(double s)
{
        double der, nJ, a, b;
        double match_pnt = phase_shift(wave, LamSq);

        double result;

        switch(qn_I) {
        case 0: {nJ = 2.; break;}
        case 1: {nJ = 1.; break;}
        case 2: {nJ = 2.; break;}
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

double omnes::kernel(double s, double sp)
{
        double result;

        if ((sign == 0) || (s <= sthPi) || (s > LamSq))
        {
                result = extrap_phase(sp) / (sp * (sp - s));
        }
        else
        {
                result = (extrap_phase(sp) - extrap_phase(s)) / (sp * (sp - s));
        }
};

std::complex<double> omnes::eval(double s)
{
        if (WG_GENERATED == false)
        {
                double weights[N_omnes], abscissas[N_omnes];
                gauleg(sthPi, LamSq, abscissas, weights, N_omnes);
                std::cout << " Generating Gaussian-Legendre Quandrature weight for Omnes function... " << "\n";

                for (int i = 0; i < N_omnes; i++)
                {
                        wgt.push_back(weights[i]);
                        abs.push_back(abscissas[i]);
                }

                WG_GENERATED = true;
        }

        std::complex<double> result;
        double sum = 0.;
        double alpha = extrap_phase(s) / M_PI;

        for (int i = 0; i < N_omnes; i++)
        {
                sum += kernel(s, abs[i]) * wgt[i];
        }

        if (sign == 0)
        {
                std::cout << "working \n";
                result = exp( (s/M_PI) * sum / pow(LamSq - s - xi*double(sign)*eps, alpha));
        }
        else
        {
                std::cout << "not working \n";
                result = 0;
        }
        return result;
};

void omnes::print_wgt()
{
        std::cout << wgt.size() << "\n";
};
