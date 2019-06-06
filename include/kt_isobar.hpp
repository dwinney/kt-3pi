// These objects include the isobars (i.e. the partial waves in each channel)
// they inherite from the amplitude class because they describe the scattering dynamics
// of a single spin and isospin-projected subchannel.
//
// Dependencies: amp.hpp, omnes.hpp, iteration.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _ISOBAR_
#define _ISOBAR_

#include <cmath>

#include "kt_iteration.hpp"
#include "iomanip"

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

  // Members related to iterating over the KT integral
  vector<iteration> iters;
  void start(double mass);
  iteration next(iteration previous);

//-----------------------------------------------------------------------------

public:
  isobar(int isospin, int spin, double mass) :
  spin_proj(spin), iso_proj(isospin), omega(isospin, spin)
  {
    start(mass);
  };

  ~isobar(){};

};
//-----------------------------------------------------------------------------



#endif
