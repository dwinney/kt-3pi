// These are methods for evaluating the KT angular (t) integral (i.e. the angular projection of cross subchannel
// isobars).
// This is seperated to allow the possibility of testing different methods of evaluating dispersion (s) integral.
//
// Dependencies: kt_iteration.hpp, decay_kinematics.hpp, aux_math.hpp, omnes.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _T_INT_
#define _T_INT_

#include "kt_iteration.hpp"
#include "decay_kinematics.hpp"
#include "aux_math.hpp"
#include "omnes.hpp"

class angular_integral
{
private:
  friend class kt_equations;
  iteration * previous;

  int N_integ = 60;

  decay_kinematics kinematics;
  double mDec = kinematics.get_decayMass();

protected:
  complex<double> k(double s); //Analytically continued k-momentum function
  complex<double> z_s(double s, double t); // Scattering angle in s-channel

  // Kinematic kernel
  // This can be extended to have isospin and helicity dependence
  complex<double> kernel(double s, double t);

  // the threshold and psuedo_threshold points for the 2 -> 2 scattering process
  double a0 = 0.5 *  (mDec * mDec - mPi * mPi);

  // Bounds of integration
  double t_minus(double s);
  double t_plus(double s);

  // Integration regions. Naming convention match Igor Danilkin's mathematica notebook for clarity
  // s0 = sthPi
  // a = pseudo_thresh
  // b = threshold
  complex<double> integ_s0_a0(double s);
  complex<double> integ_a0_a(double s);
  complex<double> integ_a_b(double s);
  complex<double> integ_b(double s);

public:
  // Default constructor
   angular_integral(decay_kinematics dec)
  : kinematics(dec) {};

  // Evaluate the inhomogenous contribution at a given energy
  complex<double> operator() (double s);
};

#endif
