// Class definitions for KT formalism, including the rescattering integral and KT amplitude
//
// Dependencies: kt_isobar.hpp
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

#include "kt_isobar.hpp"
#include "kt_equations.hpp"
#include "kt_options.hpp"

//-----------------------------------------------------------------------------
// The kt class combines isobars together.
class kt_amplitude
{
protected:
// The maximum number of iterations of the KT integral and the maximal spin_projection
kt_options options;
kt_equations kt;
void start();

// Storing all the relevant info of the decay particle
decay_kinematics kinematics;

// Storing successive iterations of the KT equations
vector<iteration> iters;

//-----------------------------------------------------------------------------
public:
kt_amplitude(kt_options ops,  decay_kinematics kine)
  : options(ops), kinematics(kine),
    kt(ops, kine)
  {};

void iterate();

// Print the nth iteration
void print_iteration(int n, int m);
//-----------------------------------------------------------------------------
};
#endif
