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
// Here we evaluate using the standard method integrating to infinity or
// the conformal mapping method used in [1409.7708].
//
// Depends on the value of use_conformal in kt_options
// ---------------------------------------------------------------------------
  class dispersion_integral
  {
  protected:
    decay_kinematics kinematics;
    kt_options options;

    const double LamOmnes = previous->omega.LamOmnes; // Conformal dispersion cutoff

    angular_integral inhom;
    const double a = kinematics.pseudo_threshold(); // Pseudo threshold, problematic point!
    const double b = kinematics.threshold();

    iteration * previous;

    //Explicitly factor out the factors of k(s) give singularities
    complex<double> Mtilde(int n, double s, int ieps);
    complex<double> reg_integrand(int n, double s, int ieps); //the regular piece of the integrand.

    // Integrate around the singularity at a
    const int N_integ = 15;
    double interval = 0.005; // Interval on either side of a to integrate

    complex<double> integrate(int n, double s, int ieps, double low, double up);

    complex<double> disperse(int n, double s, int ieps);
    complex<double> conformal_log_reg(int n, double s, int ieps);

  public:
    // Default constructor
    dispersion_integral(kt_options ops, decay_kinematics dec)
      : options(ops), inhom(ops, dec), kinematics(dec)
    { };

    complex<double> operator() (int n, double s, int ieps);

    void pass_iteration(iteration * prev);
  };

#endif
