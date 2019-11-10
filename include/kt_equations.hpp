// These are methods for evaluating the dispersion integral to give the final KT solution.
//
// Here we evaluate using the conformal mapping method used in [1409.7708].
//
// The omnes function is only evaluated in the elastic region [4 m_pi^2 , ~1] GeV and inelastic contributions
// are parameterized by a polynomial expansion in a conformal variable.
//
// The angular integral is done by explicitly by analytical continuing the momenta and integrating over
// the complex plane see "kt_ang_integral.hpp"
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _KT_EQ_
#define _KT_EQ_

#include "kt_options.hpp"
#include "kt_isobar.hpp"
#include "kt_ang_integral.hpp"
#include "kt_disp_integral.hpp"
#include "utilities.hpp"
#include "omnes.hpp"

// ---------------------------------------------------------------------------
// Each isobar (i.e. each spin and isospin projected partial wave) will have its own
// instance of the kt_equations object. This object allows the calculation of each next iterations
// of the rescattering ladder. Specific to those quantum numbers.
//
// Not callable, instead use kt_equations::iterate with some pointer to an iteration
// to output the next iteration (i.e. interpolations above and below unitarity cut and a copy of the original omnes function)
// ---------------------------------------------------------------------------

class kt_equations
{
private:
  // friend class isobar;
  friend class kt_amplitude;

  // This object holds all the relevant kinematic quantities such as masses and QN's
  decay_kinematics kinematics;

  // This object holds all user-input information
  // and a function to print the settings to command line
  kt_options options;

  // A method for evaluating the unitarity correction, i.e. the dispersion integral over the discontinuity.
  dispersion_integral dispersion;
  subtraction_polynomial polynomial;

  // In interpolalations exclude an interval around the pseudo_threshold
  double exc = 0.05;

  // calculate the jth (isobar index not spin!)
  isobar iterate_isobar(iteration * prev, int bar_index);

public:
  kt_equations(kt_options ops, decay_kinematics dec)
  : kinematics(dec), options(ops),
    dispersion(ops, dec),
    polynomial(options.use_conformal)
  {};

  // Calculate the next iteration from the previous one
  iteration iterate(iteration * prev);
};

#endif
