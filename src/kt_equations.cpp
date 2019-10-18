// Class which calculates and assembles the KT solution
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "kt_equations.hpp"

// ----------------------------------------------------------------------------
// Display function for all the settings input from the stored kt_options
void kt_equations::print_options()
{
  cout << "Using KT equations for ";
  if (kinematics.get_decayParticle() != "")
  {
    cout << kinematics.get_decayParticle() << ", ";
  }
  cout << "Mass = " << kinematics.get_decayMass() << " GeV, ";
  cout << "J^PC = " << kinematics.get_JPC() << endl;
  cout << "-> with max spin, j_max = " << options.max_spin << endl;
  cout << "-> with ";
  cout << options.max_iters << " Iterations, ";
  cout << options.max_subs << " Subtractions." << endl;
  cout << std::boolalpha << "-> with USE_CONFORMAL = " << options.use_conformal << "." << endl;
  cout << "-> and with Interpolation up to s = " << options.interp_cutoff << " GeV." << endl;
};

// ----------------------------------------------------------------------------
// Take in a pointer to a previous iteration (because all previous isobars are required too get next one)
//  and calculates above and below the unitarity cut
// according to the methods in the dispersion_integral class.
// Produces a new iteration with the new values.
iteration kt_equations::iterate(iteration * prev)
{
  vector<isobar> isobars;
  // sum over spins
  for (int i = 0; 2*i+1 <= options.max_spin; i++)
  {
    cout << " -> Calculating isobar with spin (" << 2*i+1 << "/" << options.max_spin << ")... " << endl;
    isobars.push_back(iterate_isobar(prev, i));

    // add a space in the terminal output if theres more isobars for aesthetic purposes
    if (2*(i+1)+1 <= options.max_spin)
    {
    cout << endl;
    }
  }

  iteration next(prev->N_iteration + 1, isobars);

  return next;
};

// ----------------------------------------------------------------------------
// This functions calculates the next iteration of the jth isobar.
isobar kt_equations::iterate_isobar(iteration * prev, int j)
{
  // Dont forget to pass the pointer to the dispersion intergral!!
  // disp(spin, subtraction, s, ieps) is the call to the dispersion integral
  disp.pass_iteration(prev);

  // each isobar has n subtractions that need to be solved for independently
  // here lambda = 1 and I = 1 are fixed
  isobar next(1, 2*j+1, 1, options, kinematics);

  // Calculate each 'fundamental solution' in the basis of subtractions
  for (int n = 0; n <= options.max_subs; n++)
  {
      cout << " --> Calculating subtraction (" << n << "/" << options.max_subs << ")... " << endl;

      // exclude a small interval around the singulartiy at pseudo_threshold
      // the interpolation will make up for it
      double a = kinematics.pseudo_threshold();
      double frac = ((a - exc) - (sthPi + EPS )) / (options.interp_cutoff - sthPi - EPS);
      int num = int(frac * double(interpolation::N_interp));

      vector<double> s;
      vector<complex<double>> above, below;

      // store values from threshold to a
      for (int i = 0; i < num; i ++)
      {
        double s_i = (sthPi + EPS) + double(i) * ((a - exc) - (sthPi + EPS)) / double(num - 1);
        s.push_back(s_i);

        // above unitarity cut (ieps = +1)
        complex<double> ab = poly(n, s_i, +1) + disp(j, n, s_i, +1);
        ab *= prev->isobars[j].omega(s_i, +1);
        above.push_back(ab);

        // below unitarity cut (ieps = -1)
        complex<double> be = poly(n, s_i, -1) + disp(j, n, s_i, -1);
        be *= prev->isobars[j].omega(s_i, -1);
        below.push_back(be);
      }

      // and then from a to the cutoff
      for (int i = 0; i < (interpolation::N_interp - num) ; i++)
      {
        double s_i = (a + exc) + double(i) * (options.interp_cutoff - (a + exc)) / double(interpolation::N_interp - num - 1);
        s.push_back(s_i);

        complex<double> ab = poly(n, s_i, +1) + disp(j, n, s_i, +1);
        ab *= prev->isobars[j].omega(s_i, +1);
        above.push_back(ab);

        complex<double> be = poly(n, s_i, -1) + disp(j, n, s_i, -1);
        be = prev->isobars[j].omega(s_i, -1);
        below.push_back(be);
      }

      subtraction nth(n, s, above, below);
      next.subtractions.push_back(nth);
  }

  return next;
};
