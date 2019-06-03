// These objects include the isobars (i.e. the partial waves in each channel)
// they inherite from the amplitude class because they describe the scattering dynamics
// of a single spin and isospin-projected subchannel.
//
// Dependencies: amp.hpp, omnes.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _ISOBAR_
#define _ISOBAR_

#include <iostream>
#include <complex>
#include <vector>

#include "amplitude.hpp"
#include "omnes.hpp"

#include "Math/Interpolator.h"

using std::vector;
using std::complex;
using std::cout;

// Wrapper class to better interface with ROOT's interpolation classes for both real and imaginary parts
//-----------------------------------------------------------------------------
class interpolation
{
protected:
  vector<double> s, r_fx, i_fx;

//-----------------------------------------------------------------------------
public:
  interpolation(vector<double> si, vector<complex<double>> fx)
  {
    // we need to split the complex<double> into two vector<double>
    for (int i = 0; i < 200; i++)
    {
      s.push_back(si[i]);

      double re = real(fx[i]);
      double im = imag(fx[i]);

      r_fx.push_back(re);
      i_fx.push_back(im);
    }
  };
};

// Each iteration contains the values of the current function and an interpolation
//-----------------------------------------------------------------------------
class iteration
{
protected:
  interpolation above, on, below;
//-----------------------------------------------------------------------------
public:
  // Constructor for the 0th iteration. Takes values of the bare Omnes function
  iteration(vector<double> s, vector<complex<double>> sabove, vector<complex<double>> son, vector<complex<double>> sbelow)
   : above(s, sabove), on(s, son), below(s, sbelow)
    {};

  // // Constructor for the N > 0 iteration, takes the previous one.
  // iteration(const iteration &previous){};

  ~iteration(){};

  int print(int x);
};
//-----------------------------------------------------------------------------

// Each isobar should contain the quantum numbers of the whole ampltidue
// and information on its isospin and spin projections
//-----------------------------------------------------------------------------
class isobar :  public omnes
{
protected:
  int spin_proj, iso_proj;

  int N_inter = 0;
  vector<iteration> iters;
//-----------------------------------------------------------------------------

public:
  isobar(int isospin, int spin) :
  spin_proj(spin), iso_proj(isospin), omnes(isospin, spin)
  { };
  ~isobar(){};

  void start();
  void next();
};
//-----------------------------------------------------------------------------



#endif
