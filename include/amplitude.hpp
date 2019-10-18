// Abstract class to define an amplitude to easier interact with things like the
// Dalitz object.
//
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _AMPLITUDE_
#define _AMPLITUDE_

#include "constants.hpp"

class amplitude
{
public:
  amplitude(){};

  // Store all kinematic functions and quantities for reaction
  decay_kinematics kinematics;

  // Evalute the amplitude at some energies s, t in the Dalitz region
  virtual complex<double> eval(double s, double t) = 0;
};

#endif
