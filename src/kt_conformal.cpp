// These are methods for evaluating the KT angular integral (i.e. the angular projection of cross subchannel
// isobars). This is seperated to allow the possibility of testing different methods of evaluating the integral.
//
// Here we evaluate using the conformal mapping method used in [1409.7708].
//
// The omnes function is only evaluated in the elastic region [4 m_pi^2 , ~1] GeV and inelastic contributions
// are parameterized by a polynomial expansion in a conformal variable.
// The angular integral is done by explicitly by analytical continuing the momenta and integrating over
// the complex plane.
//
// Dependencies: omnes.hpp, aux_math.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "kt_conformal.hpp"

// Analytically continued momentum k
// ---------------------------------------------------------------------------
complex<double> conformal_int::k(double s)
{
  complex<double> temp1 = sqrt(xr * (s - sthPi) / s);
  complex<double> temp2 = sqrt(xr * (threshold - s)) * sqrt(xr * (pseudo_thresh - s));

  return temp1 * temp2;
};

// COM Scattering angle in the s-channel scattering subchannel
// Complex in general because of the path of integration needed by analytic continuation
// ---------------------------------------------------------------------------
complex<double> conformal_int::z_s(double s, double t)
{
  complex<double> k_s = k(s);
  double temp1 = 2. * t + s - mDec * mDec - 3. * mPi * mPi;

  return temp1 * xr / k_s;
};

// Bounds of integration
// ---------------------------------------------------------------------------
double conformal_int::t_minus(double s)
{
  if (s < sthPi || s > (mDec + mPi) * (mDec + mPi))
  {
    cout << "t_minus: Energy outside the physical decay region! Quitting... \n";
    exit(1);
  }
  else
  {
    return (mDec * mDec + 3. * mPi * mPi - s) / 2. + real(k(s)) / 2.;
  }
};

double conformal_int::t_plus(double s)
{
  if (s < sthPi || s > (mDec + mPi) * (mDec + mPi))
  {
    cout << "t_plus: Energy outside the physical decay region! Quitting... \n";
    exit(1);
  }
  else
  {
    return (mDec * mDec + 3. * mPi * mPi - s) / 2. - real(k(s)) / 2.;
  }
};

// Kinematic kernel inside the angular integral. In general this should depend
// on the helicity and spin of the isobar.
//
// Right now it only has the omega case of Lambda = 1, j = 1, I = I
// ---------------------------------------------------------------------------
complex<double> conformal_int::kernel(double s, double t)
{
  complex<double> zs, temp1, temp2;
  zs = z_s(s, t);
  temp1 = (xr - zs * zs);
  temp2 = pow((pseudo_thresh - s), 1.5);

  return 3. * temp1 * temp2 / k(s);
};

// Integration path in the complex plane
// ---------------------------------------------------------------------------
complex<double> conformal_int::integ_s0_a0(double s)
{
    if (s < sthPi || s > a0)
    {
      cout << "integ_s0_a0: Integration out of range! Quitting... \n";
      exit(1);
    }

    double w[N_integ + 1], x[N_integ + 1];
    gauleg(t_minus(s), t_plus(s), w, x, N_integ);

    complex<double> sum = 0.;
    for (int i = 1; i < N_integ + 1; i++)
    {
      complex<double> temp = kernel(s, x[i]) * interp_above(x[i]);
      sum += w[i] * temp;
    }
    return sum;
};

complex<double> conformal_int::integ_a0_a(double s)
{
  if (s < a0|| s > pseudo_thresh)
  {
    cout << "integ_a0_a: Integration out of range! Quitting... \n";
    exit(1);
  }

  double wM[N_integ + 1], xM[N_integ + 1];
  gauleg(t_minus(s), sthPi - 1.e-5, wM, xM, N_integ);

  double wP[N_integ + 1], xP[N_integ + 1];
  gauleg(sthPi - 1.e-5, t_plus(s), wP, xP, N_integ);

  complex<double> sumP = 0., sumM = 0.;
  for(int i = 1; i < N_integ + 1; i++)
  {
    complex<double> tempM = kernel(s, xM[i]) * interp_below(xM[i]);
    complex<double> tempP = kernel(s, xP[i]) * interp_above(xP[i]);

    sumM += wM[i] * tempM;
    sumP += wP[i] * tempP;
  }

  return sumM + sumP;
};

complex<double> conformal_int::integ_a_b(double s)
{
  if (s < pseudo_thresh || s > threshold)
  {
    cout << "integ_a_b: Integration out of range! Quitting... \n";
    exit(1);
  }

  omega.set_ieps(0);
  double wM[N_integ + 1], xM[N_integ + 1];
  gauleg(t_minus(s), 2. * t_minus(threshold), wM, xM, N_integ);

  double wP[N_integ + 1], xP[N_integ + 1];
  gauleg(2. * t_plus(threshold), t_plus(s), wP, xP, N_integ);

  complex<double> sumP = 0., sumM = 0.;
  for(int i = 1; i < N_integ + 1; i++)
  {
    complex<double> tempM = kernel(s, xM[i]) * omega.eval(xM[i]);
    complex<double> tempP = kernel(s, xP[i]) * omega.eval(xP[i]);

    sumM += wM[i] * tempM;
    sumP += wP[i] * tempP;
  }

  return sumM + sumP;
};

complex<double> conformal_int::integ_b(double s)
{
  if (s < threshold)
  {
    cout << "integ_s0_a0: Integration out of range! Quitting... \n";
    exit(1);
  }

  omega.set_ieps(0);
  double w[N_integ + 1], x[N_integ + 1];
  gauleg(t_minus(s), t_plus(s), w, x, N_integ);

  complex<double> sum = 0.;
  for (int i = 1; i < N_integ + 1; i++)
  {
    complex<double> temp = kernel(s, x[i]) * omega.eval(x[i]);
    sum += w[i] * temp;
  }
  return sum;
};

// Angular integral, F_hat
// ---------------------------------------------------------------------------
complex<double> conformal_int::Fhat(double s)
{
  if (abs(s - sthPi) < 0.0001)
  {
    return 2. * interp_above(t_minus(sthPi)) * pow((pseudo_thresh - sthPi), 1.5);
  }
  else if (s > sthPi && s < a0)
  {
    return integ_s0_a0(s);
  }
  else if (s >= a0 && s <= pseudo_thresh)
  {
    return integ_a0_a(s);
  }
  else if (s > pseudo_thresh && s < threshold)
  {
    return integ_a_b(s);
  }
  else if (s >= threshold)
  {
    return integ_b(s);
  }
};

complex<double> conformal_int::Mtilde(double s)
{
  double absOmnes = abs(omega.eval(s));
  double sinDelta = sin(omega.extrap_phase(s));

  return sinDelta * Fhat(s) / absOmnes;
};

complex<double> conformal_int::solution(double s, int ieps = +1)
{
  omega.set_ieps(ieps);
  complex<double> Omnes = omega.eval(s);

  complex<double> s_integral = 
};
