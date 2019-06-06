// These are methods for evaluating the KT angular integral (i.e. the angular projection of cross subchannel
// isobars). This is seperated to allow the possibility of testing different methods of evaluating the integral.
//
// Here we evaluate using the conformal mapping method used in [1409.7708].
//
// The omnes function is only evaluated in the elastic region [4 m_pi^2 , ~1] GeV and inelastic contributions
// are parameterized by a polynomial expansion in a conformal variable.
// The angular integral is done by explicitly by analytical continuing the momenta and integrating over
// the complex plane.
//
// Dependencies: omnes.hpp, aux_math.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _CONFORMAL_
#define _CONFORMAL_

#include "aux_math.hpp"
#include "omnes.hpp"

class conformal_int
{
protected:
  int N_integ = 60;
  double LamDisp = 1.; // Dispersion cutoff.

  omnes omega;
  interpolation interp_above, interp_below;

  complex<double> k(double s); //Analytically continued k-momentum function
  complex<double> z_s(double s, double t); // Scattering angle in s-channel
  complex<double> kernel(double s, double t);

  // the threshold and psuedo_threshold points for the 2 -> 2 scattering process
  double mDec;
  double threshold = (mDec + mPi) * (mDec + mPi);
  double pseudo_thresh = (mDec - mPi) * (mDec - mPi);
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

  // The integral, combines all the above at a given s
  complex<double> Fhat(double s);
  complex<double> Mtilde(double s);
  complex<double> solution(double s, int ieps);

public:
  // Pass by reference to not make a copy of the interpolations.
  conformal_int(double mass, omnes & on, interpolation & above, interpolation & below)
  : mDec(mass), interp_above(above), interp_below(below), omega(on)
  {};
};

#endif
