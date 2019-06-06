// Internal classes used by the isobar class to incorporate 3-body effects.
// The iteration object contains the actual content of the KT equations
//
// Dependencies: amp.hpp, omnes.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _ITER_
#define _ITER_

#include "aux_math.hpp"
#include "amplitude.hpp"
#include "omnes.hpp"
#include "kt_conformal.hpp"

//-----------------------------------------------------------------------------
// Each iteration contains the values of the current function and an interpolation
// Here is where the actually "KT equations" are done
//-----------------------------------------------------------------------------
class iteration
{
protected:
  int N_iteration;

  double mDec;

  interpolation interp_above, interp_below;
  omnes omega;

//-----------------------------------------------------------------------------
public:
  // Constructor for the 0th iteration. Takes values of the bare Omnes function
  iteration(int n, double mass, omnes omeg, vector<double> s,
                                            vector<complex<double>> above,
                                            vector<complex<double>> below)
   : N_iteration(n), mDec(mass), omega(omeg), interp_above(s, above), interp_below(s, below)
    {};

  // Copy constructor needed to define vector<iterations>
  // Additionally is the constructor for the N > 0 th iteration whih takes the previous
  iteration(const iteration &previous)
  : N_iteration(previous.N_iteration),
    mDec(previous.mDec),
    omega(previous.omega),
    interp_above(previous.interp_above),
    interp_below(previous.interp_below)
  {};

  // Destructor
  ~iteration(){};

};

#endif
