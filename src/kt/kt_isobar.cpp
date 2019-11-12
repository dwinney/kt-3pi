// These objects include the isobars (i.e. the partial waves in each channel)
// they inherite from the amplitude class because they describe the scattering dynamics
// of a single spin and isospin-projected subchannel.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "kt/kt_isobar.hpp"

// ----------------------------------------------------------------------------
  // Check how many subtractions this isobar has as specified in kt_options
void isobar::check_subtractions()
{
  for (int i = 0; i < options.subIDs.size(); i++)
  {
    if (options.subIDs[i][0] == iso_proj && options.subIDs[i][1] == spin_proj && options.subIDs[i][2] == hel_proj)
    {
      n_subs = options.subIDs[i][3];
    }
  }
};

// ----------------------------------------------------------------------------
// Populate the zeroth iteration with values directly from call to Omnes function
void isobar::zeroth()
{
  // If there are previously stored iterations for some reason, clear it when start() is called.
  if (subtractions.size() != 0)
  {
    subtractions.clear();
    cout << "isobar: Clearing previously stored subtractions for j = " << spin_proj;
    cout << " and I = " << iso_proj << ". \n";
  }

  cout << "-> initializing isobar with I = " << iso_proj;
  cout << ", j = " << spin_proj;
  cout << ", and lambda = " << hel_proj;
  // cout << " with " << n_subs << " fundamental solutions...";
  cout << "... \n";

  // Loop over the basis of (n = number of subtractions) "fundamental solutions"
  subtraction_polynomial sub_poly(options.use_conformal);
  for (int n = 0; n <= n_subs; n++)
  {
    // Need to store values of the Omnes amplitude around the unitarity cut for the
    // analytic continuation of the angular integral.
    vector<double> s;
    vector< complex<double> > above_cut, below_cut;
    complex<double> ab, be;

    for (int i = 0; i < interpolation::N_interp; i++)
    {
      double s_i = sthPi + EPS + double(i) * (options.interp_cutoff - sthPi) / double(interpolation::N_interp);
      s.push_back(s_i);

      // above the unitarity cut (+ i epsilon)
      ab =  omega(s_i, +1) * sub_poly(n, s_i, +1);
      above_cut.push_back(ab);

      // below the cut (- i epsilon)
      be =  omega(s_i, -1) * sub_poly(n, s_i, -1);
      below_cut.push_back(be);
    }

    // store the zeroth iteration of each fundamental solution in the
    // subtractions vector member of the isobar object
    subtraction sub_zero(n, s, above_cut, below_cut);
    subtractions.push_back(sub_zero);
  }
};

// ----------------------------------------------------------------------------
// Evaluate the isobar partial wave at some energy s
// Sums subtractions with their coefficients
complex<double> isobar::eval_isobar(double s)
{
  complex<double> result = 0.;

  if (s > sthPi && s <= options.interp_cutoff)
  {
    // the "Zeroth" subtraction constant is constrained by the overall Normalization
    // and thus we treat this one seperately
    result += subtractions[0].interp_above(s);

    // then sum over fundamental solutions for each subtraction order and multiply by coefficient
    for (int i = 1; i <= options.max_subs; i++)
    {
      result += coefficients[i-1] * subtractions[i].interp_above(s);
    }
  }

  else
  {
    cout << "isobar: Attempting to evaluate outside of interpolation range. Quitting..." << endl;
    exit(1);
  }

  return result;
};

// ----------------------------------------------------------------------------
// Set the subtraction coefficients for fitting routines
void isobar::set_params(vector<double> par)
{
  int n_param;
  if (options.use_conformal == true){n_param = par.size();}
  else {n_param = par.size()/2;}

  if (n_param != n_subs)
  {
    cout << "isobar::set_params() number of params and coefficients dont match! Quitting... \n";
    exit(0);
  }

  // clear any existing coefficients
  coefficients.clear();
  for (int i = 0; i < n_param; i++)
  {
    complex<double> c_i;

    // filter whether using complex or real coefficients
    if (options.use_conformal == false)
    {
      c_i = par[2*i] * exp(xi * par[2*i+1]);
    }
    else
    {
      c_i = xr * par[i];
    }

    coefficients.push_back(c_i);
  }
};

// ----------------------------------------------------------------------------
// cout a list of the parameters currently saved
void isobar::print_params()
{
  for (int i = 0; i < n_subs; i++)
  {
    string varname = "c[" + std::to_string(spin_proj) + std::to_string(i + 1) + "]";

    if (options.use_conformal == false)
    {
      cout << std::left << std::setw(20) << varname + ":";
      cout << std::setw(20) << coefficients[i] << "\n";

      cout << std::left << std::setw(20) << "mod " + varname + ":";
      cout << std::setw(20) << abs(coefficients[i]) << "\n";

      cout << std::left <<  std::setw(20) << "arg " + varname + ":";
      cout << std::setw(20) << std::arg(coefficients[i]) << "\n";
    }
    else
    {
      cout << std::left << std::setw(20) << varname + ":";
      cout << std::setw(20) << abs(coefficients[i]) << "\n";
    }
  }
};
