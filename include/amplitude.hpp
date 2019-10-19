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

#include <complex>
#include "decay_kinematics.hpp"
#include "constants.hpp"

class amplitude
{
public:
  amplitude(decay_kinematics kine)
  : kinematics(kine)
  {};

  // Store all kinematic functions and quantities for reaction
  decay_kinematics kinematics;

  //
  virtual void set_params(int n, const double *par) = 0;

  // Evalute the amplitude at some energies s, t in the Dalitz region
  virtual std::complex<double> eval(double s, double t) = 0;
};

#endif
