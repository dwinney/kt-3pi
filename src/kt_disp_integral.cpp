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

// ----------------------------------------------------------------------------
// Evaluate the dispersion integral.
// isobar index = j, subtraction_ID = n
// Method used depends on use_conformal in kt_options
complex<double> dispersion_integral::operator() (int j, int n, double s, int ieps)
{
  // Error check
  if (ieps != -1 && ieps != +1)
  {
    cout << "dispersion_integral: Invalid ieps. Quitting..." << endl;
    exit(1);
  }

  return disperse(j, n, s, ieps);
};

// ---------------------------------------------------------------------------
// Disperson integral of inhomogeneity
complex<double> dispersion_integral::disperse(int j, int n, double s, int ieps)
{
  complex<double> result;

  if (options.use_conformal == true)
  {
    if (LamOmnes < a)
    {
      result = con_integrate(j, n, s, ieps, sthPi + 2. * EPS, LamOmnes) + con_sp_log(j, n, s, ieps, sthPi + EPS, LamOmnes);
      result -= cutoff_log(j, n, s, ieps);
    }
    else
    {
      result =  con_integrate(j, n, s, ieps, sthPi + 2. * EPS, a - interval) + con_sp_log(j, n, s, ieps, sthPi + EPS, a - interval);
      result += con_integrate(j, n, s, ieps, a + interval, LamOmnes) + con_sp_log(j, n, s, ieps, a + interval, LamOmnes);
      result -= cutoff_log(j, n, s, ieps);
    }
  }

  else
  {
    result =  std_integrate(j, n, s, ieps, sthPi + 2. * EPS, a - interval) + std_sp_log(j, n, s, ieps, sthPi + EPS, a - interval);
    result += std_integrate_inf(j, n, s, ieps, a + interval) + std_sp_log_inf(j, n, s, ieps, a + interval);
  }

  return result;
};

// ---------------------------------------------------------------------------
// The regular (singularity-free) piece of the integrand
complex<double> dispersion_integral::disp_function(int j, int n, double s, int ieps)
{
  complex<double> result = inhomogeneity(j, n, s);
  result /= kinematics.barrier_factor(2*j+1, 1, s);
  result *= sin(previous->isobars[j].extrap_phase(s));
  result /= std::abs(previous->isobars[j].omega(s, ieps));

  return result;
};

// ---------------------------------------------------------------------------
// INTEGRATION FUNCTIONS FOR CONFORMAL SCHEME

// ---------------------------------------------------------------------------
// integrate dispersion kernal from low to high
complex<double> dispersion_integral::con_integrate(int j, int n, double s, int ieps, double lower, double upper)
{
  double w[N_integ + 1], x[N_integ + 1];
  gauleg(lower, upper, x, w, N_integ);

  complex<double> sum = 0.;
  for (int i = 1; i < N_integ + 1; i++)
  {
    complex<double> temp = disp_function(j, n, x[i], ieps) - disp_function(j, n, s, ieps);
    temp /=  (x[i] - s - double(ieps) * xi * EPS);

    sum += w[i] * temp;
  }

  return sum / M_PI;
};

// ---------------------------------------------------------------------------
// Analytical Integral of the subtracted term which regularizes the singularity at s^prime = s
complex<double> dispersion_integral::con_sp_log(int j, int n, double s, int ieps, double lower, double upper)
{
  complex<double> result, log_term;
  result = disp_function(j, n, s, ieps);

  log_term = log(upper - s - double(ieps) * xi * EPS);
  log_term -= log(lower - s - double(ieps) * xi * EPS);

  return result * log_term / M_PI;
};

// ----------------------------------------------------------------------------
// Log term to regularize cut-off point in the conformal mapping scheme
complex<double> dispersion_integral::cutoff_log(int j, int n, double s, int ieps)
{
  complex<double> result = disp_function(j, n, LamOmnes, ieps);
  result *= log( (LamOmnes - s - xi * double(ieps) * EPS) / (LamOmnes - sthPi));

  return result / M_PI;
};

// ---------------------------------------------------------------------------
// INTEGRATION FUNCTIONS FOR STANDARD SCHEME

// ---------------------------------------------------------------------------
// integrate dispersion kernal from low to high
complex<double> dispersion_integral::std_integrate(int j, int n, double s, int ieps, double lower, double upper)
{
  double w[N_integ + 1], x[N_integ + 1];
  gauleg(lower, upper, x, w, N_integ);

  complex<double> sum = 0.;
  for (int i = 1; i < N_integ + 1; i++)
  {
    complex<double> temp = disp_function(j, n, x[i], ieps) * pow(s * xr / x[i], previous->isobars[j].n_subs);
    temp -= disp_function(j, n, s, ieps);
    temp *= s / x[i];
    temp /= (x[i] - s - double(ieps) * xi * EPS);

    sum += w[i] * temp;
  }

  return sum / M_PI;
};

// ----------------------------------------------------------------------------
// log term that regularizes the pole at s = sp above
complex<double> dispersion_integral::std_sp_log(int j, int n, double s, int ieps, double low, double high)
{
 complex<double> result = disp_function(j, n, s, ieps);

 complex<double> log_term = log(xr * (high - s - double(ieps) * xi * EPS)) - log(xr * high);
 log_term -= log(xr * (low - s - double(ieps) * xi * EPS)) - log(xr * low);

 return result * log_term / M_PI;
};

// ----------------------------------------------------------------------------
// integrate from low to infinity
complex<double> dispersion_integral::std_integrate_inf(int j, int n, double s, int ieps, double low)
{
  double w[N_integ + 1], x[N_integ + 1];
  gauleg(0., 1., x, w, N_integ);

  complex<double> sum = 0.;
  for (int i = 1; i < N_integ + 1; i++)
  {
    double sp = low + tan(M_PI * x[i] / 2.);

    complex<double> temp = disp_function(j, n, sp, ieps) * pow(xr *(s/sp), previous->isobars[j].n_subs);
    temp -= disp_function(j, n, s, ieps); // subtract away the pole at s = x[i];
    temp *= s / sp;
    temp /=  (sp - s - double(ieps) * xi * EPS);

    temp *=  (M_PI / 2.) / pow(cos(M_PI * x[i] / 2.), 2.); // jacobian
    sum += w[i] * temp;
  }

  return sum / M_PI;
};

// ----------------------------------------------------------------------------
// log term that regularizes the pole at s = sp adapted with one subtraction to be regular when integrating to infinity.
complex<double> dispersion_integral::std_sp_log_inf(int j, int n, double s, int ieps, double low)
{
 complex<double> result = disp_function(j, n, s, ieps);
 complex<double> log_term = - (log(xr * (low - s - double(ieps) * xi * EPS)) - log(xr * low));

 return result * log_term / M_PI;
};

// ----------------------------------------------------------------------------
// Calculate the sum_rule value of the second subtraction coefficient
complex<double> dispersion_integral::sum_rule(iteration * prev)
{
    pass_iteration(prev);

    double x[N_integ + 1], w[N_integ + 1];
    gauleg(sthPi + EPS, a - interval, x, w, N_integ);
    complex<double> sum_1 = 0.;
    for (int i = 1; i < N_integ + 1; i++)
    {
      complex<double> temp = disp_function(0, 0, x[i], +1);
      temp /= (x[i] * x[i]);
      sum_1 += w[i] * temp;
    }

    double w2[N_integ + 1], x2[N_integ + 1];
    gauleg(0., 1., x2, w2, N_integ);

    complex<double> sum_2 = 0.;
    for (int i = 1; i < N_integ + 1; i++)
    {
      double sp = a + interval + tan(M_PI * x2[i] / 2.);

      complex<double> temp = disp_function(0, 0, sp, 1);
      temp /= sp * sp;
      temp *=  (M_PI / 2.) / pow(cos(M_PI * x2[i] / 2.), 2.); // jacobian
      sum_2 += w2[i] * temp;
    }

    return (sum_1 + sum_2) / M_PI;
};

//----------------------------------------------------------------------------
// UTILITY FUNCTIONS

// ----------------------------------------------------------------------------
// utility function that updates the stored pointer to the current iteration
void dispersion_integral::pass_iteration(iteration * prev)
{
    previous = NULL;

    previous = prev;
    inhomogeneity.pass_iteration(prev);
};
