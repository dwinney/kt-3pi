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
// Dependencies: kt_iteration, kt_ang_integral, aux_math, omnes
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _KT_EQ_
#define _KT_EQ_

#include "kt_iteration.hpp"
#include "kt_options.hpp"
#include "kt_ang_integral.hpp"
#include "kt_disp_integral.hpp"
#include "aux_math.hpp"
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
  // This object holds all the relevant kinematic quantities such as masses and QN's
  decay_kinematics kinematics;
  kt_options options;

  // A method for evaluating the unitarity correction, i.e. the dispersion integral over the discontinuity.
  dispersion_integral disp;
  subtraction_polynomial poly;

public:
  // TODO: KT equations depend on spin projection and helicity in general
  kt_equations(decay_kinematics dec, kt_options ops)
  : kinematics(dec), options(ops), disp(ops, dec),
    poly(options.use_conformal)
  {
    cout << "Using KT equations for ";
    if (dec.get_decayParticle() != "")
    {
      cout << dec.get_decayParticle() << ", ";
    }
    cout << "Mass = " << dec.get_decayMass() << " GeV, ";
    cout << "J^PC = " << dec.get_JPC() << endl;
    cout << "-> with ";
    cout << options.max_iters << " Iterations, ";
    cout << options.max_subs << " Subtractions." << endl;
    cout << std::boolalpha << "-> with USE_CONFORMAL = " << options.use_conformal << "." << endl;
   };

  // Calculate the next iteration from the previous one
  iteration iterate(iteration * prev);
};

#endif
