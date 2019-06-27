// Class which calculates and assembles the KT solution
//
// Dependencies: kt_iteration, kt_ang_integral, aux_math, omnes
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "kt_equations.hpp"

// ----------------------------------------------------------------------------
// Take in a pointer to a previous iteration and calculates above and below the unitarity cut
// according to the methods in the dispersion_integral class.
// Produces a new iteration with the new values.
iteration kt_equations::iterate(iteration * prev)
{
  disp.pass_iteration(prev);

  vector<subtraction> subs;

  for (int n = 0; n < prev->max_subs + 1; n++)
  {
      vector<double> s;
      vector<complex<double>> above, below;
      for (int i = 0; i < interpolation::N_interp; i++)
      {
        double s_i = sthPi + EPS + double(i) * (omnes::LamOmnes - sthPi - EPS) / double(interpolation::N_interp);
        s.push_back(s_i);

        complex<double> ab = prev->omega(s_i, +1) * (poly(n, s_i) + disp(n, s_i, +1));
        above.push_back(ab);

        complex<double> be = prev->omega(s_i, -1) * (poly(n, s_i) + disp(n, s_i, -1));
        below.push_back(be);
      }

      subtraction nth(n, s, above, below);
      subs.push_back(nth);
  }

  //Save the values to a new iteration. Add one to the internal iteration number.
  iteration next(prev->N_iteration + 1, prev->max_subs, prev->omega, subs);

  return next;
};
