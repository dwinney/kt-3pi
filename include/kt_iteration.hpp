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

#include "kt_isobar.hpp"

//-----------------------------------------------------------------------------
// Container class; each iteration holds a vector containing the current interpolated
// values of the isobars. All isobars are required to calculate next iteration
//-----------------------------------------------------------------------------
class iteration
{
public:
  // TODO: may be extended to take in 3 vectors for isospin = 0, 1, 2
  iteration(int n, vector<isobar> isos)
   : N_iteration(n), isobars(isos)
    {};

  iteration(const iteration &previous)
  : N_iteration(previous.N_iteration),
    isobars(previous.isobars)
  {};

  // Destructor
  ~iteration(){};

  const int N_iteration; // iteration ID
  vector<isobar> isobars;
};
#endif
