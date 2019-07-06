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
// Mtilde has the factors of k(s) in the inhomogeneitys explicitly factored out
complex<double> dispersion_integral::Mtilde(int n, double s, int ieps)
{
  return inhom(n, s) * pow(inhom.k(s), 3.); // Power of three because I = 1
};

// ---------------------------------------------------------------------------
// The regular (singularity-free) piece of the integrand
// TODO: Here there would be a power of s to the number of subtractions in general
complex<double> dispersion_integral::reg_integrand(int n, double s, int ieps)
{
  return sin(previous->omega.extrap_phase(s)) * Mtilde(n, s, ieps) / std::abs(previous->omega(s, ieps));;
};

// ---------------------------------------------------------------------------
// Integrals on either side of the singularity at a.
// These are discontinuity-free and we can integrate with no worries
complex<double> dispersion_integral::integ_sthPi_a(int n, double s, int ieps)
{
  double w[N_integ/2 + 1], x[N_integ/2 + 1];
  gauleg(sthPi + EPS, a - interval, x, w, N_integ/2);

  complex<double> sum = 0.;
  for (int i = 1; i < N_integ/2 + 1; i++)
  {
    complex<double> temp = reg_integrand(n, x[i], ieps) / pow(inhom.k(x[i]), 3.);
    temp -= reg_integrand(n, s, ieps) / pow(inhom.k(s), 3.); // subtract the pole at s = x[i];
    temp /=  (x[i] - s - double(ieps) * xi * EPS);

    // Here we add powers of s^prime from the subtractions in the standard scheme
    // they surpress the high energy contributions and help the integral converge.
    if (options.use_conformal == false)
    {
      temp /= pow(x[i] * xr, (options.max_subs + 1));
    }

    sum += w[i] * temp;
  }

  sum += (reg_integrand(n, s, ieps) /  pow(inhom.k(s), 3.))
          * log( (a - interval - s - double(ieps) * xi * EPS)
                / (sthPi + EPS - s - double(ieps) * xi * EPS) );

  // Matching powers of s to cancel mass dimension of subtractions
  if (options.use_conformal == false)
  {
    sum *= pow(s * xr, (options.max_subs + 1)) / M_PI;
  }

  return sum;
};

// Integrate from a to some upper limit.
// If in the conformal scheme, upper_limit = LamOmnes
// else its off at infinity (or sufficiently large i.e. 600.)
complex<double> dispersion_integral::integ_a_Lam(int n, double s, int ieps)
{
  double upper_limit;
  if (options.use_conformal == true)
  {
    upper_limit = LamOmnes;
  }
  else
  {
    upper_limit = 600.;
  }

  double w[N_integ/2 + 1], x[N_integ/2 + 1];
  gauleg(a + interval, LamOmnes, x, w, N_integ/2);

  complex<double> sum = 0.;
  for (int i = 1; i < N_integ/2 + 1; i++)
  {
    complex<double> temp = reg_integrand(n, x[i], ieps) / pow(inhom.k(x[i]), 3.);
    temp -= reg_integrand(n, s, ieps) / pow(inhom.k(s), 3.); // subtract away the pole at s = x[i];
    temp /=  (x[i] - s - double(ieps) * xi * EPS);

    // Here we add powers of s^prime from the subtractions in the standard scheme
    if (options.use_conformal == false)
    {
      temp /= pow(x[i] * xr, (options.max_subs + 1));
    }

    sum += w[i] * temp;
  }

  sum += (reg_integrand(n, s, ieps) /  pow(inhom.k(s), 3.))
          * log( (LamOmnes - s - double(ieps) * xi * EPS)
                / (a + interval - s - double(ieps) * xi * EPS) );

  // Matching powers of s to cancel mass dimension of subtractions
  if (options.use_conformal == false)
  {
    sum *= pow(s * xr, (options.max_subs + 1)) / M_PI;
  }

  return sum;
};

// ---------------------------------------------------------------------------
// In the conformal mapping scheme we only integrate from [sthPi, LamOmnes],
// so if the psuedo threshold is outside of the interval then we dont have to worry
// and we just integrate regularly through the entire region
complex<double> dispersion_integral::integ_sthPi_Lam(int n, double s, int ieps)
{
  double w[N_integ + 1], x[N_integ + 1];
  gauleg(sthPi + EPS, LamOmnes, x, w, N_integ);

  complex<double> sum = 0.;
  for (int i = 1; i < N_integ + 1; i++)
  {
    complex<double> temp = reg_integrand(n, x[i], ieps) / pow(inhom.k(x[i]), 3.);
    temp -= reg_integrand(n, s, ieps) / pow(inhom.k(s), 3.); // subtract the pole at s = x[i];
    temp /=  (x[i] - s - double(ieps) * xi * EPS);

    sum += w[i] * temp;
  }

  sum += (reg_integrand(n, s, ieps) /  pow(inhom.k(s), 3.))
          * log( (LamOmnes - s - double(ieps) * xi * EPS)
                / (a + interval - s - double(ieps) * xi * EPS) );
  return sum;
};


// ---------------------------------------------------------------------------
// Disperson integral of inhomogeneity
complex<double> dispersion_integral::disperse(int n, double s, int ieps)
{
  complex<double> result;
  if ((options.use_conformal == true) && (a < LamOmnes))
  {
    result = integ_sthPi_Lam(n, s, ieps);
  }
  else
  {
    result =  integ_sthPi_a(n, s, ieps) + integ_a_Lam(n, s, ieps);
  }

  return result;
};

// ----------------------------------------------------------------------------
// Log term to regularize cut-off point in the conformal mapping scheme
complex<double> dispersion_integral::log_reg(int n, double s, int ieps)
{
  complex<double> mtilde = sin(previous->omega.extrap_phase(LamOmnes)) * inhom(n, LamOmnes) / abs(previous->omega(LamOmnes, ieps));
  complex<double> endpoint = mtilde / pow(a * xr - LamOmnes, 1.5);
  endpoint *= log( (LamOmnes - s - xi * double(ieps) * EPS) / (LamOmnes - sthPi));

  return endpoint;
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
  return disperse(n, s, ieps) - log_reg(n, s, ieps);
  }

  else
  {
    return disperse(n, s, ieps);
  }
};
