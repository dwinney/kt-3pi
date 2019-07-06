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
#include "kt_options.hpp"

//-----------------------------------------------------------------------------
// The kt class combines isobars together.
class kt_amplitude
{
protected:
// The maximum number of iterations of the KT integral and the maximal spin_projection
kt_options options;

// Storing all the relevant info of the decay particle
decay_kinematics kinematics;

// Here the isobars are stored. Only one dimensional vector for now (omega case), need to make more when isospin is involved.
vector<isobar> isobars;

//-----------------------------------------------------------------------------
public:
kt_amplitude(kt_options ops,  decay_kinematics kine)
  : options(ops)
  {};

//-----------------------------------------------------------------------------
};
#endif
