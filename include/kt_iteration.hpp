// Internal classes used by the isobar class to incorporate 3-body effects.
// The iteration object contains the actual content of the KT equations
//
// Dependencies: amp.hpp, omnes.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _ITER_
#define _ITER_

#include "aux_math.hpp"
#include "amplitude.hpp"
#include "omnes.hpp"

#include "Math/Interpolator.h"

//-----------------------------------------------------------------------------
// Wrapper class to better interface with ROOT's interpolation class
// for both real and imaginary parts
//-----------------------------------------------------------------------------
class interpolation
{
protected:
  vector<double> s, r_fx, i_fx;
  ROOT::Math::Interpolator r_inter, i_inter;

public:
  // The number of interpolation points used  can be changed here
  const static int N_interp = 200;

  interpolation(vector<double> x, vector<complex<double>> fx)
  : s(x), r_fx(vec_real(fx)), i_fx(vec_imag(fx)),
    r_inter(x, vec_real(fx), ROOT::Math::Interpolation::kCSPLINE),
    i_inter(x, vec_imag(fx), ROOT::Math::Interpolation::kCSPLINE)
  { };

  // Copy Constructor (needed for vector<iterations>)
  // Additionally, ROOT's Interpolator object cannot be copied so this will create new
  // Interpolator with the same data.
  interpolation(const interpolation &old_interp)
  : s(old_interp.s), r_fx(old_interp.r_fx), i_fx(old_interp.i_fx),
    r_inter(old_interp.s, old_interp.r_fx, ROOT::Math::Interpolation::kCSPLINE),
    i_inter(old_interp.s, old_interp.i_fx, ROOT::Math::Interpolation::kCSPLINE)
    {};

  // The output, at given s, outputs the real and imaginary parts from the interpolations
  complex<double> operator ()(double s);
};

//-----------------------------------------------------------------------------
// Each iteration contains the values of the current function and an interpolation
// Here is where the actually "KT equations" are done
//-----------------------------------------------------------------------------
class iteration
{
protected:
  int N_iteration, N_integ = 60;
  interpolation interp_above, interp_below;
  omnes omega;

  double mDec;

  complex<double> k(double s); //Analytically continued k-momentum function
  complex<double> z_s(double s, double t); // Scattering angle in s-channel
  complex<double> kernel(double s, double t);

  // the threshold and psuedo_threshold points for the 2 -> 2 scattering process
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

//-----------------------------------------------------------------------------
public:
  // Constructor for the 0th iteration. Takes values of the bare Omnes function
  iteration(int n, double mass, omnes omeg, vector<double> s,
                                            vector<complex<double>> above,
                                            vector<complex<double>> below)
   : N_iteration(n), mDec(mass), omega(omeg), interp_above(s, above), interp_below(s, below)
    {};

  // Copy constructor needed to define vector<iterations>
  // Additionally is the constructor for the N > 0 th iteration whih takes the previous
  iteration(const iteration &previous)
  : N_iteration(previous.N_iteration),
    mDec(previous.mDec),
    omega(previous.omega),
    interp_above(previous.interp_above),
    interp_below(previous.interp_below)
  {};

  // Destructor
  ~iteration(){};


};

#endif
