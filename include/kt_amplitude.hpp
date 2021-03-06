// Amplitude class that assembles all isobars togethers.
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

double normalization = 1.;

void start();

//-----------------------------------------------------------------------------
public:
// Constructor
kt_amplitude(kt_options ops,  decay_kinematics kine)
  : amplitude(kine),
    options(ops), kt(ops, kine)
  {};

// Storing successive iterations of the KT equations
vector<iteration> iters;

void iterate();

// Normalize to some experimental value
void normalize(double gamma_exp);

// Evaluate the total amplitude at some energys s and t
complex<double> eval(double s, double t);
void set_params(int n, const double *par)
{};

// Print the nth iteration
void plot_iteration(int n, int j, int m);
void plot_isobar(int n);

double error_func(double s, double t)
{
  return 1.;
};

//-----------------------------------------------------------------------------
};
#endif
