// These are methods for evaluating the KT dispersion (s) integral.
// This is seperated from the angular integral, and the total KT equations which combine everything
// because of the delicate nature of handleing the singularities in the integrand near pseudo threshold.
//
// Dependencies: kt_iteration.hpp, decay_kinematics.hpp, aux_math.hpp, omnes.hpp
//               angular_integral.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _DISP_
#define _DISP_

#include "kt_iteration.hpp"
#include "kt_ang_integral.hpp"
#include "decay_kinematics.hpp"
#include "aux_math.hpp"
#include "omnes.hpp"

// ---------------------------------------------------------------------------
// Callable by operator() outputs the unitary correction to the omnes amplitude
// associated with rescattering at some real energy s, and either above or below
// unitarity cut ieps = +1 or -1.
//
// Here we evaluate using the conformal mapping method used in [1409.7708].
//
// The omnes function is only evaluated in the elastic region [4 m_pi^2 , ~1] GeV and inelastic contributions
// are parameterized by a polynomial expansion in a conformal variable.
// ---------------------------------------------------------------------------

class dispersion_integral
{
protected:
  decay_kinematics kinematics;
  const double LamOmnes = previous->omega.LamOmnes; // Dispersion cutoff

  angular_integral inhom;
  const double a = inhom.a; // Pseudo threshold, problematic point!

  iteration * previous;

  //Explicitly factor out the factors of k(s) give singularities
  complex<double> Mtilde(double s, int ieps);
  complex<double> reg_integrand(double s, int ieps); //the regular piece of the integrand.

  // Integral from s = [sthPi, LamOmnes] split into three pieces:
  // 1: from [sthPi, a - interval]
  // 2: from [a - interval, a + interval]
  // 3: from [a + interval, LamOmnes]

  const int N_integ = 60;
  double interval = 0.005; // Interval on either side of a to integrate

  // End cap integrals
  complex<double> integ_sthPi_a(double s, int ieps); // 1
  complex<double> integ_a_Lam(double s, int ieps); // 3


  complex<double> disp_inhom(double s, int ieps);

  // Log term subtracted to regularize the cut off at LamOmnes
  complex<double> log_reg(double s, int ieps);

public:
  // Default constructor
  dispersion_integral(decay_kinematics dec)
    : inhom(dec), kinematics(dec)
  {};

  complex<double> operator() (double s, int ieps);

  void pass_iteration(iteration * prev);
};

#endif
