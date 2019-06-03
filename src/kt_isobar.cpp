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

#include "kt_isobar.hpp"

void isobar::start()
{
  vector<double> s;
  vector< complex<double> > above_cut, below_cut, on_cut;

  for (int i = 0; i < 200; i++)
  {
    double s_i = sthPi + double(i) * (elastic_cutoff - sthPi) / 200.;
    s.push_back(s_i);

    set_ieps(1); // above the unitarity cut (+ i epsilon)
    above_cut.push_back(omnes::eval(s_i));

    set_ieps(0); // on the cut
    on_cut.push_back(omnes::eval(s_i));

    set_ieps(-1); // below the cout
    below_cut.push_back(omnes::eval(s_i));
  }

  iteration first(s, above_cut, on_cut, below_cut);
  iters.push_back(first);
};
