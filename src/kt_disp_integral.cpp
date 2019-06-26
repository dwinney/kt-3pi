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
complex<double> dispersion_integral::Mtilde(double s, int ieps)
{
  return inhom(s) * pow(inhom.k(s), 3.); // Power of three because I = 1
};

// ---------------------------------------------------------------------------
// The regular (singularity-free) piece of the integrand
// TODO: Here there would be a power of s to the number of subtractions in general
complex<double> dispersion_integral::reg_integrand(double s, int ieps)
{
  return sin(previous->omega.extrap_phase(s)) * Mtilde(s, ieps) / std::abs(previous->omega(s, ieps));;
};

// ---------------------------------------------------------------------------
// Integrals on either side of the singularity interval.
// These are discontinuity free and we can integrate with no worries
complex<double> dispersion_integral::integ_sthPi_a(double s, int ieps)
{
  double w[N_integ/2 + 1], x[N_integ/2 + 1];
  gauleg(sthPi + EPS, a - interval, x, w, N_integ/2);

  complex<double> sum = 0.;
  for (int i = 1; i < N_integ/2 + 1; i++)
  {
    complex<double> temp = reg_integrand(x[i], ieps) / pow(inhom.k(x[i]), 3.);
    temp -= reg_integrand(s, ieps) / pow(inhom.k(s), 3.); // subtract the pole at s = x[i];
    temp /=  (x[i] - s - double(ieps) * xi * EPS);

    sum += w[i] * temp;
  }

  sum += (reg_integrand(s, ieps) /  pow(inhom.k(s), 3.))
          * log( (a - interval - s - double(ieps) * xi * EPS)
                / (sthPi + EPS - s - double(ieps) * xi * EPS) );

  return sum;
};

complex<double> dispersion_integral::integ_a_Lam(double s, int ieps)
{
  double w[N_integ/2 + 1], x[N_integ/2 + 1];
  gauleg(a + interval, LamOmnes, x, w, N_integ/2);

  complex<double> sum = 0.;
  for (int i = 1; i < N_integ/2 + 1; i++)
  {
    complex<double> temp = reg_integrand(x[i], ieps) / pow(inhom.k(x[i]), 3.);
    temp -= reg_integrand(s, ieps) / pow(inhom.k(s), 3.); // subtract the pole at s = x[i];
    temp /=  (x[i] - s - double(ieps) * xi * EPS);

    sum += w[i] * temp;
  }

  sum += (reg_integrand(s, ieps) /  pow(inhom.k(s), 3.))
          * log( (LamOmnes - s - double(ieps) * xi * EPS)
                / (a + interval - s - double(ieps) * xi * EPS) );

  return sum;
};


// ---------------------------------------------------------------------------
// Disperson integral of inhomogeneity
complex<double> dispersion_integral::disp_inhom(double s, int ieps)
{
  complex<double> result =  integ_sthPi_a(s, ieps) + integ_a_Lam(s, ieps);

  if (s > LamOmnes || s <= sthPi)
  {
    cout << "dispersion_integral: Trying to evaluate outside of Omnes range! Quitting..." << endl;
    exit(1);
  }

  return result;
};

// ----------------------------------------------------------------------------
// Log term to regularize cut-off point
complex<double> dispersion_integral::log_reg(double s, int ieps)
{
  complex<double> mtilde = sin(previous->omega.extrap_phase(LamOmnes)) * inhom(LamOmnes) / abs(previous->omega(LamOmnes, ieps));
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
// Log term to regularize cut-off point
complex<double> dispersion_integral::operator() (double s, int ieps)
{
  if (ieps != -1 && ieps != +1)
  {
    cout << "dispersion_integral: Invalid ieps. Quitting..." << endl;
    exit(1);
  }

  return disp_inhom(s, ieps) - log_reg(s, ieps);
};
