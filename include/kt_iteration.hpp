// Internal classes used by the isobar class to incorporate 3-body effects.
// The iteration object contains the actual content of the KT equations
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
// Each iteration contains the values of the current function and an interpolation
// of the KT solution above and below the unitarity cut in the complex plane.
//
// For in between threshold and pseudo-threshold on the real axis, we also save a copy of
// the omnes function (TODO: Save an interpolation instead of the whole thing)
//-----------------------------------------------------------------------------
class iteration
{
protected:


//-----------------------------------------------------------------------------
public:
  // Default constructor with no parameters passed
  // Needed for kt_equations to be defined only by decay_kinematics
  iteration(){};

  // Constructor for the 0th iteration. Takes values of the bare Omnes function
  iteration(int n, omnes omeg,  vector<double> s,
                                vector<complex<double>> above,
                                vector<complex<double>> below)
   : N_iteration(n), omega(omeg), interp_above(s, above), interp_below(s, below)
    {};

  // Copy constructor needed to define vector<iterations>
  // Additionally is the constructor for the N > 0 th iteration whih takes the previous
  iteration(const iteration &previous)
  : N_iteration(previous.N_iteration),
    omega(previous.omega),
    interp_above(previous.interp_above),
    interp_below(previous.interp_below)
  {};

  // Destructor
  ~iteration(){};

  int N_iteration;
  interpolation interp_above, interp_below;
  omnes omega;


};

#endif
