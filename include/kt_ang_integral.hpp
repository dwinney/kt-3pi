// These are methods for evaluating the KT angular (t) integral (i.e. the angular projection of cross subchannel
// isobars).
// This is seperated to allow the possibility of testing different methods of evaluating dispersion (s) integral.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _T_INT_
#define _T_INT_

#include "kt_sub_polynomial.hpp"
#include "kt_options.hpp"
#include "kt_iteration.hpp"
#include "decay_kinematics.hpp"
#include "omnes.hpp"

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TError.h>

#include <iomanip>
#include <fstream>

using std::setw;

// ---------------------------------------------------------------------------
// Internal class used by dispersion_integral. Callable through operator().
// Outputs the angular average of the previous iteration as a fixed value of
// (real) energy s.
//
// The path of integration is distorted in the complex plane and integration is done
// in regions split according to the analytic continuation of the momentum k(s).
// ---------------------------------------------------------------------------

class angular_integral
{
private:
  friend class dispersion_integral;

  iteration * previous; // Pointer to current iteration

  const int N_angular;

  decay_kinematics kinematics;
  kt_options options;
  subtraction_polynomial sub_poly;

  double mDec = kinematics.get_decayMass();
  double a0 = 0.5 *  (mDec * mDec - mPi * mPi);
  double a = kinematics.pseudo_threshold();
  double b = kinematics.threshold();

protected:
  complex<double> k(complex<double> s); //Analytically continued k-momentum function
  complex<double> sine_sqr(double s, complex<double> t); // Scattering angle in s-channel

  // Kinematic kernel
  // This can be extended to have isospin and helicity dependence
  complex<double> kernel(int j, int lam, int jp, int lamp, double s, complex<double> t);

  // Bounds of integration
  complex<double> t_minus(double s);
  complex<double> t_plus(double s);

  // Integration regions. Naming convention match Igor Danilkin's mathematica notebook for clarity
  // a = pseudo_thresh
  // b = threshold
  complex<double> integ_sthPi_a0(int j, int n, double s);
  complex<double> integ_a0_a(int j, int n, double s);
  complex<double> integ_a_b(int j, int n, double s);
  complex<double> integ_b(int j, int n, double s);

public:
  // Default constructor
   angular_integral(kt_options ops, decay_kinematics dec)
  : options(ops), kinematics(dec), sub_poly(ops.use_conformal),
    N_angular(ops.N_ang)
   {};

  // Evaluate the inhomogenous contribution at a given energy
  complex<double> operator() (int j, int n, double s);

  void pass_iteration(iteration * prev);
};

#endif
