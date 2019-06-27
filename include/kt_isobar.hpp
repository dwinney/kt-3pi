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

#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include <TStyle.h>
#include "TError.h"

#include "decay_kinematics.hpp"
#include "kt_iteration.hpp"
#include "kt_equations.hpp"
#include "iomanip"

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
  int spin_proj, iso_proj;
  omnes omega;

  // Vector storing each iteration of the KT equation
  int num_subtractions;
  vector<double> coefficients;
  vector<iteration> iters;

  // Start() to populate the 0th vector entry with interpolations of the base omnes function
  void start();

  // KT equations object
  kt_equations kt;
//-----------------------------------------------------------------------------

public:
  isobar(int isospin, int spin, int subtracts, decay_kinematics & dec) :
  spin_proj(spin), iso_proj(isospin), num_subtractions(subtracts),
  kinematics(dec), omega(isospin, spin), kt(dec)
  { };

  decay_kinematics kinematics;

  void iterate(int n);

  // Print the nth iteration
  void print(int n, int m);

  // These functions are to interface with dalitz_fit
  void set_params(int n_params, const double *par);

  // Evaluate the isobar in one channel or the total amplitude
  complex<double> subtracted_isobar(double s);
  complex<double> eval(double s, double t);
};
//-----------------------------------------------------------------------------



#endif
