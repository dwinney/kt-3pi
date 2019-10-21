// General purpose classes for Dalitz plot generation.
// for fitting to polynomial expansion see param_fit.hpp
//
// Dependencies: constants.hpp aux_math.hpp, ROOT
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _DALITZ_
#define _DALITZ_

#include "amplitude.hpp"
#include "constants.hpp"
#include "utilities.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using std::vector;
using std::setw;
using std::cout;
using std::endl;

//-----------------------------------------------------------------------------
// Define a Dalitz plot object with some other amplitude object containing a model.
//
// Initiate a dalitz plot object with:
// model_amp my_model;
// dalitz my_dalitz(&my_model);
//-----------------------------------------------------------------------------

class dalitz
{
protected:
  amplitude * amp;

//-----------------------------------------------------------------------------
// For integration
  int n = 60; // Number of integration points
  int N_int(){
    return n;
  };

  bool S_WG_GENERATED, T_WG_GENERATED;
  vector<double> s_wgt, s_abs;
  vector<vector<double>> t_wgt, t_abs; //"two-dimensional" vectors
  void generate_s_weights();
  void generate_t_weights(vector<double> s);
  void generate_weights();

  // Calculate the area of the dalitz region
  double dalitz_area();
  double offset = 0.00001;
//-----------------------------------------------------------------------------
public:
  // Default Constructor
  dalitz(amplitude * my_amp)
  : amp(my_amp)
  {};

//--------------------------------------------------------------------------
  // Double differential cross section
  double d2Gamma(double s, double t);
  double decay_width();

  // Center of the dalitz region
  double s_c = amp->kinematics.s_c();
  double t_c = amp->kinematics.t_c();
//--------------------------------------------------------------------------
  // Print out a txt file with the Dalitz plot and plot with ROOT
  void plot(string options = "");

  void set_integration_points(int n);
};

//-----------------------------------------------------------------------------
#endif
