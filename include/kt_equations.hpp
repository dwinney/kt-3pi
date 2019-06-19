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
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _KT_EQ_
#define _KT_EQ_

#include "kt_iteration.hpp"
#include "kt_ang_integral.hpp"
#include "aux_math.hpp"
#include "omnes.hpp"

// TODO: figure out how to only have to generate Gaussian weights once
class kt_equations
{
private:
  int N_integ = 60;

  // This object holds all the relevant kinematic quantities such as masses and quantu n
  decay_kinematics kinematics;
  const double pthresh = kinematics.pseudo_threshold();
  const double mDec = kinematics.get_decayMass();

  // Evaluation of the angular integral in the complex plane
  inhomogeneity inhom;

protected:
  // Pointer to the previous iteration used to evaluate the next one
  iteration * previous;

  double LamDisp = 1.; // Dispersion cutoff.
  complex<double> disp_inhom(double s, int ieps); // Disperson integral of inhomogeneity
  complex<double> log_reg(double s, int ieps); //Log term to regularize elastic cuttoff
  complex<double> solution(double s, int ieps);

public:
  // TODO: KT equations depend on spin projection and helicity in general
  kt_equations(decay_kinematics dec)
  : kinematics(dec), inhom(dec)
  {
    cout << "Using KT equations for ";
    if (dec.get_decayParticle() != "")
    {
      cout << dec.get_decayParticle() << ", ";
    }
    cout << "Mass = " << dec.get_decayMass() << " GeV, ";
    cout << "J^PC = " << dec.get_JPC() << endl;
   };

  // Calculate the next iteration from the previous one
  iteration iterate(iteration * prev);
};

#endif
