// Amplitude class that assembles all isobars togethers.
// Inherits from the amplitude class and can therefore be used with dalitz and dalitz_fit
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _KT_
#define _KT_

#include <iostream>
#include <complex>
#include <vector>

#include "amplitude.hpp"
#include "kt_isobar.hpp"
#include "kt_equations.hpp"
#include "kt_options.hpp"
#include "dalitz.hpp"

//-----------------------------------------------------------------------------
// The kt class combines isobars together.
// Requires a kt_options object which contains all relevant information on
// subtractions and spin sums, and a decay_kinematics for kinematic quantities
//-----------------------------------------------------------------------------

class kt_amplitude : public amplitude
{
private:
// The maximum number of iterations of the KT integral and the maximal spin_projection
kt_options options;
kt_equations kt;
void print_options(); // print command line summary of all options recieved

double normalization = 1.; // overall normalization

// initialize with bare omnes solutions
void start();

//-----------------------------------------------------------------------------
public:
// Constructor
kt_amplitude(kt_options ops,  decay_kinematics kine)
  : amplitude(kine),
    options(ops), kt(ops, kine)
  {
    print_options();
    start();
  };

// Calculate and store successive iterations of the KT equations
vector<iteration> iters;
void iterate();

// Normalize to some experimental value
void normalize(double gamma_exp);

// Set and print subtraction coeffs for fitting
void set_params(int n, const double *par);
void print_params();
void sum_rule();

// Evaluate the total amplitude at some energys s and t
complex<double> eval(double s, double t);


// Printing to file
void plot_fundamentalSolutions(int j);
void plot_isobar(int n, string file = "");

//-----------------------------------------------------------------------------
};
#endif
