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

    // const double LamOmnes = previous->omega.LamOmnes; // Conformal dispersion cutoff

    angular_integral inhom;
    const double LamOmnes = 1.;
    const double a = kinematics.pseudo_threshold(); // Pseudo threshold, problematic point!
    const double b = kinematics.threshold();

    iteration * previous;

    // function being dispersed, ihomogeneity * sin(delta) / |Omega|
    complex<double> disp_function(int j, int n, double s, int ieps);

    // Integrate around the singularity at a
    const int N_integ = 30;
    double interval = 0.0007; // Interval on either side of a to integrate

    // Integration functions for the conformal scheme
    complex<double> con_integrate(int j, int n, double s, int ieps, double low, double high); // integrate dispersion kernel from low to high
    complex<double> con_sp_log(int j, int n, double s, int ieps, double low, double high); // log term from subtracted pole at s = sprime
    complex<double> cutoff_log(int j, int n, double s, int ieps); // log term to regularie the cutoff at LamOmnes

    // integration for the standard scheme
    // these need to take into account the number of subtractions in the integration kernel
    // where the above do not
    complex<double> std_integrate(int j, int n, double s, int ieps, double low, double high);
    complex<double> std_integrate_inf(int j, int n, double s, int ieps, double low); // integrate from low up to inf
    complex<double> std_sp_log(int j, int n, double s, int ieps, double low, double high);
    complex<double> std_sp_log_inf(int j, int n, double s, int ieps, double low);

    // Evaluate the dispersion integral at some s
    complex<double> disperse(int j, int n, double s, int ieps);

    // utility to print a dat file and plot of the inhomogneity
    void angular_test(int j, int n);

  public:
    // Default constructor
    dispersion_integral(kt_options ops, decay_kinematics dec)
      : options(ops), inhom(ops, dec), kinematics(dec)
    { };

    complex<double> operator() (int j, int n, double s, int ieps);

    void pass_iteration(iteration * prev);

    // complex<double> sum_rule(iteration * prev);
  };

#endif
