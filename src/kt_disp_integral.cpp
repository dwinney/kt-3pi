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


complex<double> dispersion_integral::integrate(int n, double s, int ieps, double low, double up)
{

  complex<double> log_piece;
  if (options.use_conformal == true)
  {
    log_piece = (reg_integrand(n, s, ieps) /  pow(inhom.k(s), 3.))
            * log( (up - s - double(ieps) * xi * EPS)
                  / (low - s - double(ieps) * xi * EPS) );
  }
  else
  {
    //TODO THIS LOG PIECE REQUIRES GENERALIZATION
    log_piece = pow(s * xr, options.max_subs) * (reg_integrand(n, s, ieps) /  pow(inhom.k(s), 3.));
    log_piece *= log( (up - s - double(ieps) * xi * EPS) / (low - s - double(ieps) * xi * EPS))
                  - (options.max_subs + 1) * log(up / low);
  }

  double w[N_integ + 1], x[N_integ + 1];
  gauleg(low, up, x, w, N_integ);

  complex<double> sum = 0.;
  for (int i = 1; i < N_integ + 1; i++)
  {
    complex<double> temp = reg_integrand(n, x[i], ieps) / pow(inhom.k(x[i]), 3.);
    temp -= reg_integrand(n, s, ieps) / pow(inhom.k(s), 3.); // subtract away the pole at s = x[i];
    temp /=  (x[i] - s - double(ieps) * xi * EPS);

    // Here we add powers of s^prime from the subtractions in the standard scheme
    if (options.use_conformal == false)
    {
      temp *= pow(s * xr, (options.max_subs + 1)) / pow(x[i] * xr, (options.max_subs + 1));
    }

    sum += w[i] * temp;
  }

  sum += log_piece;

  return sum / M_PI;
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
    result = integrate(n, s, ieps, sthPi + EPS, LamOmnes);
    }
    else
    {
      result =  integrate(n, s, ieps, sthPi + EPS, a - interval);
      result += integrate(n, s, ieps, a + interval, LamOmnes);
    }
  }

  else
  {
    result =  integrate(n, s, ieps, sthPi + EPS, a - interval);
    result += integrate(n, s, ieps, a + interval, 600.);
  }

  return result;
};

// ----------------------------------------------------------------------------
// Log term to regularize cut-off point in the conformal mapping scheme
complex<double> dispersion_integral::conformal_log_reg(int n, double s, int ieps)
{
  complex<double> mtilde = sin(previous->omega.extrap_phase(LamOmnes)) * inhom(n, LamOmnes) / abs(previous->omega(LamOmnes, ieps));
  complex<double> endpoint = mtilde / pow(a * xr - LamOmnes, 1.5);
  endpoint *= log( (LamOmnes - s - xi * double(ieps) * EPS) / (LamOmnes - sthPi));

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
  return disperse(n, s, ieps) - conformal_log_reg(n, s, ieps);
  }

  else
  {
    return disperse(n, s, ieps);
  }
};
