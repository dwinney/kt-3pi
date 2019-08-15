// These objects include the isobars (i.e. the partial waves in each channel)
// they inherite from the amplitude class because they describe the scattering dynamics
// of a single spin and isospin-projected subchannel.
//
// Dependencies: decay_kinematics.hpp, iteration.hpp, kt_equation.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _ISOBAR_
#define _ISOBAR_

#include <cmath>
#include <fstream>
#include <iomanip>

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TError.h>

#include "kt_options.hpp"
#include "kt_sub_polynomial.hpp"
#include "decay_kinematics.hpp"
#include "omnes.hpp"

using std::setw;

//-----------------------------------------------------------------------------
// The subtraction object contains an interpolation of the current iteration of the
// KT equations for a given order of subtraction polynomial above and below
// the unitarity cut
//-----------------------------------------------------------------------------
class subtraction
{
public:
  subtraction(int n, vector<double> s,
                      vector<complex<double>> above,
                      vector<complex<double>> below)
    : N_subtraction(n), interp_above(s, above), interp_below(s, below)
    {};

  subtraction(const subtraction &previous)
  : N_subtraction(previous.N_subtraction), interp_above(previous.interp_above),
                                           interp_below(previous.interp_below)
  {};

  const int N_subtraction; // subtraction ID
  interpolation interp_above, interp_below;
};

//-----------------------------------------------------------------------------
// Each isobar should contain the quantum numbers of the whole ampltidue
// and information on its isospin and spin projections.
//
// The isobar class is the actual object that is called and combined by the KT amplitude.
//-----------------------------------------------------------------------------

class isobar : public omnes
{
protected:

  friend class kt_equations;
  friend class angular_integral;
  friend class dispersion_integral;

  // Vector storing each iteration of the KT equation
  kt_options options;

  // Subtraction coefficients (complex in general)
  double normalization = 1.;
  vector<complex<double>> coefficients;

//-----------------------------------------------------------------------------
public:
  isobar(int isospin, int spin, int helicity, kt_options opti, decay_kinematics dec) :
  spin_proj(spin), iso_proj(isospin), hel_proj(helicity),
  options(opti),
  kinematics(dec), omnes(isospin, spin, opti.use_conformal)
  {};

  isobar(const isobar & prev)
  : spin_proj(prev.spin_proj), iso_proj(prev.iso_proj), hel_proj(prev.hel_proj),
    options(prev.options),
    kinematics(prev.kinematics), omnes(prev),
    subtractions(prev.subtractions)
  {};

  ~isobar(){};

  int spin_proj, iso_proj, hel_proj;
  decay_kinematics kinematics;

  // Basis of functions for the subtracted solutions
  vector<subtraction> subtractions;

  // store the bare omnes function with the current quantum numbers I & j
  void zeroth();

  // void print();

  // These functions are to interface with dalitz_fit
  void set_params(int n_params, const double *par);
  void print_params();

  void normalize(double gamma_exp);
  void sum_rule();

  double error_func(double s, double t)
  {
    return 1.;
  };

  // Evaluate the isobar in one channel or the total amplitude
  complex<double> subtracted_isobar(double s);
};
//-----------------------------------------------------------------------------



#endif
