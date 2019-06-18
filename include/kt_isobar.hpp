// These objects include the isobars (i.e. the partial waves in each channel)
// they inherite from the amplitude class because they describe the scattering dynamics
// of a single spin and isospin-projected subchannel.
//
// Dependencies: decay_kinematics.hpp, omnes.hpp, iteration.hpp
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
  decay_kinematics kinematics;

  int spin_proj, iso_proj;
  omnes omega;

  // Vector storing each iteration of the KT equation
  vector<iteration> iters;

  // Start() to populate the 0th vector entry with interpolations of the base omnes function
  void start();

  // KT equations object
  kt_equations kt;
//-----------------------------------------------------------------------------

public:
  isobar(int isospin, int spin, decay_kinematics & dec) :
  spin_proj(spin), iso_proj(isospin), kinematics(dec), omega(isospin, spin), kt(dec)
  { };

  ~isobar(){};

  void iterate(int n);

  //Print the nth iteration
  void print(int n);
};
//-----------------------------------------------------------------------------



#endif
