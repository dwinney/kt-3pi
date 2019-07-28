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
isobar kt_equations::iterate_isobar(iteration * prev, int j)
{
  disp.pass_iteration(prev);

  // Save and interpolate only in the Dalitz region as needed for fits
  double interp_cutoff;
  if (omnes::LamOmnes > kinematics.pseudo_threshold())
  {
    interp_cutoff = omnes::LamOmnes;
  }
  else
  {
    interp_cutoff = kinematics.pseudo_threshold();
  }

  // each isobar has n subtractions that need to be solved for independently
  isobar next(1, 2*j+1, 1, options, kinematics);

  for (int n = 0; n < options.max_subs + 1; n++)
  {
      cout << " ---> Calculating subtraction (" << n << "/" << options.max_subs << ")... " << endl;

      // exclude a small interval around the singulartiy at pseudo_threshold
      double a = kinematics.pseudo_threshold();
      double frac = ((a - exc) - (sthPi + EPS )) / (interp_cutoff - sthPi + EPS);
      int num = int(frac * double(interpolation::N_interp));

      vector<double> s;
      vector<complex<double>> above, below;

      // store values from threshold to a
      for (int i = 0; i < num; i ++)
      {
        double s_i = (sthPi + EPS) + double(i) * ((a - exc) - (sthPi + EPS)) / double(num - 1);
        s.push_back(s_i);

        complex<double> ab = prev->isobars[j].omega(s_i, +1) * (poly(n, s_i, +1) + disp(j, n, s_i, +1));
        above.push_back(ab);

        complex<double> be = prev->isobars[j].omega(s_i, -1) * (poly(n, s_i, -1) + disp(j, n, s_i, -1));
        below.push_back(be);
      }

      // then from a to the cutoff
      for (int i = 0; i < (interpolation::N_interp - num) ; i++)
      {
        double s_i = (a + exc) + double(i) * (interp_cutoff - (a + exc)) / double(interpolation::N_interp - num - 1);
        s.push_back(s_i);

        complex<double> ab = prev->isobars[j].omega(s_i, +1) * (poly(n, s_i, +1) + disp(j, n, s_i, +1));
        above.push_back(ab);

        complex<double> be = prev->isobars[j].omega(s_i, -1) * (poly(n, s_i, -1) + disp(j, n, s_i, -1));
        below.push_back(be);
      }

      subtraction nth(n, s, above, below);
      next.subtractions.push_back(nth);
  }

  return next;
};

iteration kt_equations::iterate(iteration * prev)
{
  vector<isobar> isobars;
  for (int i = 0; 2*i+1 <= options.max_spin; i++)
  {
    cout << " --> Calculating isobar  with spin (" << 2*i+1 << "/" << options.max_spin << ")... " << endl;
    isobars.push_back(iterate_isobar(prev, i));
  }

  iteration next(prev->N_iteration + 1, isobars);

  return next;
};
