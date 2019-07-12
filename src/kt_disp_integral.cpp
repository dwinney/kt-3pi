// These are methods for evaluating the KT dispersion (s) integral.
// This is seperated from the angular integral, and the total KT equations which combine everything
// because of the delicate nature of handleing the singularities in the integrand near pseudo threshold.
//
// Dependencies: kt_iteration.hpp, decay_kinematics.hpp, aux_math.hpp, omnes.hpp
//               angular_integral.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "kt_disp_integral.hpp"

// ---------------------------------------------------------------------------
// The regular (singularity-free) piece of the integrand
// TODO: Here there would be a power of s to the number of subtractions in general
complex<double> dispersion_integral::reg_integrand(int n, double s, int ieps)
{
  complex<double> mtilde = inhom(n, s); // Power of three because I = 1
  return sin(previous->omega.extrap_phase(s)) * mtilde / std::abs(previous->omega(s, ieps));;
};

complex<double> dispersion_integral::integrate(int n, double s, int ieps, double low, double up)
{
  double w[N_integ + 1], x[N_integ + 1];
  gauleg(low, up, x, w, N_integ);

  complex<double> sum = 0.;
  for (int i = 1; i < N_integ + 1; i++)
  {
    complex<double> temp = reg_integrand(n, x[i], ieps) / pow(inhom.k(x[i]), 3.);

    // Here we add powers of s^prime from the subtractions in the standard scheme
    if (options.use_conformal == false)
    {
      temp *= pow(s * xr, (options.max_subs + 1)) / pow(x[i] * xr, (options.max_subs + 1));
    }

    // Subtract off F(s) to remove pole at x[i] = s
    temp -= reg_integrand(n, s, ieps) / pow(inhom.k(s), 3.); // subtract away the pole at s = x[i];
    temp /=  (x[i] - s - double(ieps) * xi * EPS);

    sum += w[i] * temp;
  }

  return sum / M_PI;
};

// ---------------------------------------------------------------------------
// Analytical Integral of the subtracted term which regularizes the singularity at s^prime = s
complex<double> dispersion_integral::sp_log(int n, double s, int ieps, double low, double up)
{
  complex<double> result, log_term;
  result = reg_integrand(n, s, ieps);
  result /=  pow(inhom.k(s), 3.);

  log_term = log(up - s - double(ieps) * xi * EPS);
  log_term -= log(low - s - double(ieps) * xi * EPS);

  return result * log_term / M_PI;
};

// ---------------------------------------------------------------------------
// Disperson integral of inhomogeneity
complex<double> dispersion_integral::disperse(int n, double s, int ieps)
{
  complex<double> result;
  // If pseudothreshold is outside of the cutoff integrate without worry
  // no discontinuities
  if (options.use_conformal == true)
  {
    if (LamOmnes < a)
    {
    result = integrate(n, s, ieps, sthPi + EPS, LamOmnes) + sp_log(n, s, ieps, sthPi + EPS, LamOmnes);
    }
    else
    {
      result =  integrate(n, s, ieps, sthPi + EPS, a - interval) + sp_log(n, s, ieps, sthPi + EPS, a - interval);
      result += integrate(n, s, ieps, a + interval, LamOmnes) + sp_log(n, s, ieps, a + interval, LamOmnes);
    }
  }

  else
  {
    complex<double> result1 =  integrate(n, s, ieps, sthPi + EPS, a - interval) + sp_log(n, s, ieps, sthPi + EPS, a - interval);
    complex<double> result2 = integrate(n, s, ieps, a + interval, cutoff);
    complex<double> result3 = sp_log(n, s, ieps, a + interval, cutoff);
      // cout << std::left << setw(15) << s - a << setw(30) << result2 << setw(30) << result3 << endl;

    result = result1 + result2 + result3;
  }

  return result;
};

// ----------------------------------------------------------------------------
// Log term to regularize cut-off point in the conformal mapping scheme
complex<double> dispersion_integral::cutoff_log(int n, double s, double cutoff, int ieps)
{
  complex<double> mtilde = sin(previous->omega.extrap_phase(cutoff)) * inhom(n, cutoff) / abs(previous->omega(cutoff, ieps));
  complex<double> endpoint = mtilde / pow(a * xr - cutoff, 1.5);
  endpoint *= log( (cutoff - s - xi * double(ieps) * EPS) / (cutoff - sthPi));

  return endpoint / M_PI;
};

void dispersion_integral::pass_iteration(iteration * prev)
{
    previous = NULL;

    previous = prev;
    inhom.pass_iteration(prev);
};

// ----------------------------------------------------------------------------
// Evaluate the dispersion integral.
// Method used depends on use_conformal in kt_options
complex<double> dispersion_integral::operator() (int n, double s, int ieps)
{
  if (ieps != -1 && ieps != +1)
  {
    cout << "dispersion_integral: Invalid ieps. Quitting..." << endl;
    exit(1);
  }
  if (options.use_conformal == true)
  {
    return disperse(n, s, ieps) - cutoff_log(n, s, LamOmnes, ieps);
  }
  else
  {
    return disperse(n, s, ieps);
  }
};
