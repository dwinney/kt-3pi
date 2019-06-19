// These are methods for evaluating the dispersion integral to give the final KT solution.
//
// Here we evaluate using the conformal mapping method used in [1409.7708].
//
// The omnes function is only evaluated in the elastic region [4 m_pi^2 , ~1] GeV and inelastic contributions
// are parameterized by a polynomial expansion in a conformal variable.
//
// The angular integral is done by explicitly by analytical continuing the momenta and integrating over
// the complex plane see "kt_ang_integral.hpp"
//
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "kt_equations.hpp"

// ---------------------------------------------------------------------------
// Disperson integral of inhomogeneity
complex<double> kt_equations::disp_inhom(double s, int ieps)
{
  previous->omega.set_ieps(ieps);

  //integrate from thresh to p-thresh
  double w1[N_integ + 1], x1[N_integ / 2 + 1];
  gauleg(sthPi + EPS, pthresh - EPS, w1, x1, N_integ / 2);
  //then from p-thresh to disperive cutoff
  double w2[N_integ/2 + 1], x2[N_integ / 2 + 1];
  gauleg(pthresh + EPS, LamDisp, w2, x2, N_integ / 2);

  complex<double> sum1, sum2;
  for (int i = 1; i < N_integ / 2 + 1; i++)
  {
    sum1 += w1[i] * sin(previous->omega.extrap_phase(x1[i])) * inhom(x1[i]) / abs(previous->omega.eval(x1[i]));
    sum1 /= (pow(pthresh - x1[i] + xi * EPS, 1.5) * (x1[i] - s - xi * double(ieps) * EPS));


    sum2 += w2[i] * sin(previous->omega.extrap_phase(x2[i])) * inhom(x2[i]) / abs(previous->omega.eval(x2[i]));
    sum2 /= (pow(pthresh - x2[i] + xi * EPS, 1.5) * (x2[i] - s - xi * double(ieps) * EPS));
  }

  return (sum1 + sum2) / M_PI;
};

complex<double> kt_equations::log_reg(double s, int ieps)
{
  previous->omega.set_ieps(ieps);

  complex<double> mtilde = sin(previous->omega.extrap_phase(LamDisp)) * inhom(LamDisp) / abs(previous->omega.eval(LamDisp));
  complex<double> endpoint = mtilde / pow(pthresh * xr - LamDisp, 1.5);
  endpoint *= log( (LamDisp - s - xi * double(ieps) * EPS) / (LamDisp - sthPi));

  return endpoint;
};

complex<double> kt_equations::solution(double s, int ieps)
{
  previous->omega.set_ieps(ieps);

  return previous->omega.eval(s) * (xr + disp_inhom(s,ieps) - log_reg(s, ieps));
};

// ----------------------------------------------------------------------------
// Take in a pointer to a previous iteration and calculates above and below the unitarity cut.
iteration kt_equations::iterate(iteration * prev)
{
  previous = prev;
  inhom.previous = prev;

  vector<double> s;
  vector<complex<double>> above, below;
  for (int i = 0; i < interpolation::N_interp; i++)
  {
    double s_i = sthPi + EPS + double(i) * (elastic_cutoff - sthPi - EPS) / double(interpolation::N_interp);
    s.push_back(s_i);

    complex<double> ab = solution(s_i, 1);
    above.push_back(ab);

    complex<double> be = solution(s_i, -1);
    below.push_back(be);
  }

  iteration next(previous->N_iteration + 1, previous->omega, s, above, below);

  return next;
};
