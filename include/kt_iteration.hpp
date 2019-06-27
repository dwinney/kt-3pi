// Container classes used by the isobar class to store successive iterations of the kt equations.
//
// Dependencies: aux_math.hpp, omnes.hpp, decay_kinematics.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _ITER_
#define _ITER_

#include "aux_math.hpp"
#include "decay_kinematics.hpp"
#include "omnes.hpp"

//-----------------------------------------------------------------------------
// The subtraction object contains an interpolation of the current iteration of the
// KT equations for a given order of subtraction polynomial
//-----------------------------------------------------------------------------
class subtraction
{
public:
  subtraction(int n, vector<double> s,
                      vector<complex<double>> above,
                      vector<complex<double>> below)
    : N_subtraction(n), interp_above(s, above), interp_below(s, below)
    {};

  subtraction(const subtraction &previous)
  : N_subtraction(previous.N_subtraction), interp_above(previous.interp_above),
                                           interp_below(previous.interp_below)
  {};

  const int N_subtraction; // subtraction ID
  interpolation interp_above, interp_below;
};

//-----------------------------------------------------------------------------
// Each iteration contains a vector containing the current interpolated values of
// the fundamental solutions for each subtration of the KT equations.
// Additionally it has a copy of the original omnes object for ease of access.
//-----------------------------------------------------------------------------
class iteration
{
public:
  iteration(int n, int m, omnes omeg,  vector<subtraction> subs)
   : N_iteration(n), max_subs(m), omega(omeg), subtractions(subs)
    {
      if (subs.size() > m + 1)
      {
        cout << "iteration: Number of subtractions greater than max in iteration decleration. Quittng..." << endl;
        exit(1);
      }
    };

  iteration(const iteration &previous)
  : N_iteration(previous.N_iteration),
    max_subs(previous.max_subs),
    omega(previous.omega),
    subtractions(previous.subtractions)
  {};

  // Destructor
  ~iteration(){};

  const int N_iteration; // iteration ID

  const int max_subs; //
  vector<subtraction> subtractions;
  omnes omega;
};

#endif
