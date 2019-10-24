// These objects include the isobars (i.e. the partial waves in each channel)
// they inherite from the amplitude class because they describe the scattering dynamics
// of a single spin and isospin-projected subchannel.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "kt_isobar.hpp"

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

// // ----------------------------------------------------------------------------
// // Set second subtraction coefficient to its sum rule values
// void isobar::sum_rule()
// {
//   //check there is two Subtractions
//   if (iters.back().subtractions.size() != 2)
//   {
//     cout << "isobar: need two subtractions to calcualte sum rule values. Quitting... \n";
//     exit(1);
//   }
//
//   cout << "Calculating sum rule value for subtraction constant..." << endl;
//
//   complex<double> b = kt.disp.sum_rule(&iters.back());
//   if (coefficients.size() != 0)
//   {
//     coefficients.clear();
//   };
//   coefficients.push_back(b);
//
//   cout << "Sum rule constant = " << abs(b) <<  " exp(" << arg(b) << "*I)" << endl;
// };

// ----------------------------------------------------------------------------
// Print the current stored subtraction coefficients in command line
void isobar::print_params()
{
    cout << "Printing Subtraction Coefficients... \n";
    cout << "---------------------------------------- \n";
    for (int i = 0; i < options.max_subs; i++)
    {
      cout << std::left <<  setw(20) << "|a_" + std::to_string(i) + "|:" <<
      setw(20) << std::abs(coefficients[i]) << "\n";
      cout << std::left <<  setw(20) << "arg a_" + std::to_string(i) + ":" <<
      setw(20) << std::arg(coefficients[i]) << "\n";
    }
    cout << "---------------------------------------- \n";
    cout << "\n";
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

//
// // ----------------------------------------------------------------------------
// // Set the subtraction coefficients for fitting routines
// void isobar::set_params(int n_params, const double *par)
// {
//   if ((n_params - 1) / 2 != options.max_subs)
//   {
//     cout << "isobar: Number of parameters and subtrations don't match! Quitting..." << endl;
//     exit(1);
//   }
//
//   // clear existing coefficients
//   coefficients.clear();
//   normalization = par[0];
//   for (int i = 0; i < (n_params - 1)/2; i++)
//   {
//     complex<double> c_i = par[2*i+1] * exp(xi * par[2*i+2]);
//     coefficients.push_back(c_i);
//   }
//
//   // check again for good measure
//   if (coefficients.size() != (n_params-1) / 2)
//   {
//   cout << "isobar: Number of parameters and subtrations don't match! Quitting..." << endl;
//   exit(1);
//   }
// };
