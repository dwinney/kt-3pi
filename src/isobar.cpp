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

#include "isobar.hpp"

// Populate the zeroth iteration with values directly from call to Omnes function
// ---------------------------------------------------------------------------
void isobar::start(double mass)
{
  vector<double> s;

  // Need to store values of the Omens amplitude around the unitarity cut for the
  // analytic continuation of the angular integral.
  vector< complex<double> > above_cut, below_cut;
  complex<double> ab, be;

  for (int i = 0; i < interpolation::N_interp; i++)
  {
    double s_i = sthPi + 0.0001 + double(i) * (elastic_cutoff - sthPi - 0.0001) / double(interpolation::N_interp);
    s.push_back(s_i);

    omega.set_ieps(1); // above the unitarity cut (+ i epsilon)
    ab = omega.eval(s_i);
    above_cut.push_back(ab);

    omega.set_ieps(-1); // below the cut (- i epsilon)
    be = omega.eval(s_i);
    below_cut.push_back(be);
  }

  // If there are previously stored iterations for some reason, clear it when start() is called.
  if (iters.size() != 0)
  {
    iters.clear();
    cout << "isobar: Clearing previously stored Iterations for j = " << spin_proj;
    cout << " and I = " << iso_proj << ". \n";
  }

  iteration zeroth(0, mass, omega, s, above_cut, below_cut);
  iters.push_back(zeroth);

};
