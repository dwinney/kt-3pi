// Functions for Omnes function procedures to solve unitarity integral equations.
//
// Dependencies: pipi
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#include "omnes.hpp"

constexpr double omnes::LamOmnes;

//-----------------------------------------------------------------------------
// Set number of Gaussian Quandrature points in evaluating Omnes function
void omnes::set_N_omnes(int i)
{
        N_omnes = i;
        cout << "omnes: Number of integration points set to " << N_omnes << "...\n";
};
//-----------------------------------------------------------------------------
// Smoothly extrapolated Phase-shift shift matched at Lambda_phase^2
// Based off Danilkin et al, phi -> 3pi analysis
// http://cgl.soic.indiana.edu/jpac/w3pi.php
double omnes::extrap_phase(double s)
{
  double der, nJ, a, b;
  double match_pnt = phase_shift(wave, LamSq);

  double result;

  switch(pipi_qn_I) {
  case 0: {nJ = 2.; break;}
  case 1: {nJ = 1.; break;}
  case 2: {nJ = 2.; break;}
  };


  if  (s <= LamSq)
  {
          result = phase_shift(wave, s);
  }
  else
  {
          der = phase_shift(wave, LamSq + hD) - phase_shift(wave, LamSq - hD);
          der /= 2.*hD;
          a = 3. * pow(nJ*M_PI - match_pnt,2) / (2. * LamSq * der);
          b = 1. + 3. * (nJ*M_PI - match_pnt) / (2. * LamSq * der);
          result = nJ * M_PI - a / (b + pow(s/LamSq, 1.5));
  }

  return result;
};

//-----------------------------------------------------------------------------
complex<double> omnes::kernel(complex<double> s, double sp, int ieps)
{
  complex<double> result;
  if (ieps == 0)
  {
          result = extrap_phase(sp) / (sp * (sp - s));
  }
  else
  {
    double s_real = real(s);
    if ((s_real <= sthPi) || (s_real > LamSq))
    {
          result = extrap_phase(sp) / (sp * (sp - s_real - double(ieps) * xr * eps));
    }
    else
    {
          result = (extrap_phase(sp) - extrap_phase(s_real)) / (sp * (sp - s_real - double(ieps) * xr * eps));
    }
  }

  return result;
};

//-----------------------------------------------------------------------------
void omnes::check_weights()
{
  if (WG_GENERATED == false)
  {

    wgt.clear(); abs.clear();

    double weights[N_omnes + 1], abscissas[N_omnes + 1];
    gauleg(0., 1., abscissas, weights, N_omnes + 1);
    // cout << "omnes: Generating Gaussian-Legendre Quandrature weights and abscissas... \n";

    for (int i = 1; i < N_omnes + 1; i++)
    {
            wgt.push_back(weights[i]);
            abs.push_back(abscissas[i]);
    }

    WG_GENERATED = true;
  }
};

//-----------------------------------------------------------------------------
complex<double> omnes::omega_0(complex<double> s, int ieps)
{
  check_weights();

  complex<double> sum = 0.;
  for (int i = 0; i < N_omnes; i++)
  {
      if (use_conformal == true)
      {
        double sp = (1. - abs[i]) * (sthPi + EPS) + abs[i] * LamOmnes;
        complex<double> temp = kernel(s, sp, ieps);
        temp *= LamOmnes -  (sthPi + EPS);

        sum += temp * wgt[i];
      }
      else
      {
        double sp = (sthPi + EPS) + tan( M_PI * abs[i] / 2.);
        complex<double> temp =  kernel(s, sp, ieps);
        temp *= (M_PI / 2.) / pow(cos(M_PI * abs[i] / 2.), 2.);

        sum += temp * wgt[i];
      }
  }

  return exp((s /M_PI) * sum);
};

//-----------------------------------------------------------------------------
// Evaluation of the omnes function with the cutoff at the elastic region
complex<double> omnes::omega_el(complex<double> s, int ieps)
{
  // if evaluating on the cut or with an imaginary s
  if (ieps == 0 && ieps == 0)
  {
          complex<double> alpha_0 = extrap_phase(LamOmnes) / M_PI;
          return omega_0(s, ieps) / pow(LamOmnes - s, alpha_0);
  }
  else if (ieps == -1 || ieps == 1)
  {
    double s_real = real(s);

    complex<double> alpha = extrap_phase(s_real) / M_PI;
    complex<double> alpha_0 = extrap_phase(LamOmnes) / M_PI;

    // unphysical values
    if (s_real <= sthPi || s_real > LamOmnes)
    {
      return omega_0(s_real, ieps) / pow(LamOmnes - s_real - double(ieps) * xi * EPS, alpha_0);
    }

    // physical values
    else
    {
      complex<double> result;
      result = pow(s0, alpha) * pow(std::abs(LamOmnes * (s_real - sthPi)), - alpha);
      result *= omega_0(s, ieps);
      result *= exp(xi * double(ieps) * extrap_phase(s_real));
      if ( std::abs(s_real - LamOmnes) > 0.001)
      {
              result *= pow( std::abs(LamOmnes - s_real), alpha - alpha_0);
      }
      return result;
    }
  }
  else
  {
    cout << "omnes: Cannot evaluate omnes function here... Quitting" << endl;
    exit(1);
  }
};
//-----------------------------------------------------------------------------
complex<double> omnes::omega_prime(complex<double> s, int ieps)
{
  if (ieps == 0)
  {
    return omega_0(s, ieps);
  }
  else if (ieps == -1 || ieps == 1)
  {
    double s_real = real(s);

    // unphysical values
    if (s_real <= sthPi || s_real > LamOmnes)
    {
      return omega_0(s_real, ieps);
    }
    else
    {
      complex<double> result;
      result = omega_0(s_real, ieps);
      result *= exp(log(s0 / (s_real - s0)) * extrap_phase(s_real) / M_PI);
      result *= exp(double(ieps) * xi * extrap_phase(s_real));

      return result;
    }
  }
};

//-----------------------------------------------------------------------------
complex<double> omnes::operator ()(complex<double> s, int ieps)
{
  if (use_conformal == true)
  {
    return omega_el(s, ieps);
  }
  else
  {
    return omega_prime(s, ieps);
  }
};
