// Functions for Omnes function procedures.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#ifndef _OMNES_
#define _OMNES_

#include "pipi.hpp"
#include "utilities.hpp"

using std::endl;

//-----------------------------------------------------------------------------
// Define an omnes object with a given partial wave and isospin:
// omnes my_omnes(isospin, spin);
//
// Evaluate the Omnes function (as a complex number) at some s:
// my_omnes.eval(s);
// -----------------------------------------------------------------------------

class omnes : public pipi
{
//-----------------------------------------------------------------------------
private:
  // Matching energy above which phase shifts are extrapolated
  static constexpr double Lambda_phase = 1.3;
  double LamSq = Lambda_phase * Lambda_phase;

  // Gaussian weights and abscissas for integration
  int N_omnes = 60;
  bool WG_GENERATED = false;
  vector<double> wgt, x;
  void check_weights();

  // Evaluation parameters
  int wave; // partial wave spin
  bool use_conformal = false; // whether to evaluate with conformal cutoff method

  // Internal functions to evaluate dispersion integral
  complex<double> kernel(complex<double> s, double sp, int ieps);
  complex<double> omega_0(complex<double> s, int ieps);
  complex<double> omega_std(complex<double> s, int ieps);
  complex<double> omega_con(complex<double> s, int ieps);

//-----------------------------------------------------------------------------
public:
// Default constructor
  omnes() : pipi()
  {
    check_weights();
  };

  // Parameterized constructor with quantum numbers
  omnes(int i, int j, bool conformal = false)
  : pipi(i), wave(j), use_conformal(conformal)
  {
    check_weights();
  };

  // Copy constructor
  omnes(const omnes &previous):
    pipi(previous.pipi::pipi_qn_I), wave(previous.wave),
    wgt(previous.wgt), x(previous.x), WG_GENERATED(previous.WG_GENERATED),
    use_conformal(previous.use_conformal)
  {};

  static constexpr double LamOmnes = 1.0;

  // Utility functions
  void set_N_omnes(int i);

  // Smoothly extrapolated phase up to infinity
  double extrap_phase(double s);

  // User-end function to evaluate the omnes function
  // at some complex s if ieps == 0
  // or real s with ieps == +- 1
  complex<double> omega(complex<double> s, int ieps);
};
//-----------------------------------------------------------------------------
#endif
