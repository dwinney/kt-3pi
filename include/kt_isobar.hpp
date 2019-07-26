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
#include "decay_kinematics.hpp"
#include "kt_iteration.hpp"
#include "kt_equations.hpp"

using std::setw;

//-----------------------------------------------------------------------------
// Each isobar should contain the quantum numbers of the whole ampltidue
// and information on its isospin and spin projections.
//
// The isobar class is the actual object that is called and combined by the KT amplitude.
//-----------------------------------------------------------------------------

class isobar
{
protected:
  int spin_proj, iso_proj, helicity_proj;
  omnes omega;

  // Vector storing each iteration of the KT equation
  kt_options options;
  vector<complex<double>> coefficients;
  vector<iteration> iters;

  // Start() to populate the 0th vector entry with interpolations of the base omnes function
  void start();

  // KT equations object
  kt_equations kt;

//-----------------------------------------------------------------------------
public:
  isobar(int isospin, int spin, int helicity, kt_options opti, decay_kinematics dec) :
  spin_proj(spin), iso_proj(isospin), helicity_proj(helicity),
  options(opti),
  kinematics(dec), omega(isospin, spin, opti.use_conformal), kt(dec, opti)
  {
    coefficients.push_back(1.);
  };

  decay_kinematics kinematics;

  void iterate();

  // Print the nth iteration
  void print_iteration(int n, int m);
  void print();

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
  complex<double> eval(double s, double t);
};
//-----------------------------------------------------------------------------



#endif
