// These are methods for evaluating the KT angular (t) integral (i.e. the angular projection of cross subchannel
// isobars). This is seperated to allow the possibility of testing different methods of evaluating dispersion (s) integral.
//
// Dependencies: omnes.hpp, aux_math.hpp, deay_kinematics.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "kt_ang_integral.hpp"

// ---------------------------------------------------------------------------
// Analytically continued momentum k
complex<double> inhomogeneity::k(double s)
{
  complex<double> temp1 = sqrt(xr * (s - sthPi) / s);
  complex<double> temp2 = sqrt(xr * (kinematics.threshold() - s)) *
                          sqrt(xr * (kinematics.pseudo_threshold() - s));

  return temp1 * temp2;
};

// ---------------------------------------------------------------------------
// COM Scattering angle in the s-channel scattering subchannel
// Complex in general because of the path of integration needed by analytic continuation
complex<double> inhomogeneity::z_s(double s, double t)
{
  complex<double> k_s = k(s);
  double temp1 = 2. * t + s - mDec * mDec - 3. * mPi * mPi;

  return temp1 * xr / k_s;
};

// ---------------------------------------------------------------------------
// Bounds of integration
double inhomogeneity::t_minus(double s)
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

double inhomogeneity::t_plus(double s)
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

// ---------------------------------------------------------------------------
// Kinematic kernel inside the angular integral. In general this should depend
// on the helicity and spin of the isobar.
//
// Right now it only has the omega case of Lambda = 1, j = 1, I = I
complex<double> inhomogeneity::kernel(double s, double t)
{
  complex<double> zs, temp1, temp2;
  zs = z_s(s, t);
  temp1 = (xr - zs * zs);
  temp2 = pow((kinematics.pseudo_threshold() - s), 1.5);

  return 3. * temp1 * temp2 / k(s);
};

// ---------------------------------------------------------------------------
// Integration path in the complex plane
complex<double> inhomogeneity::integ_s0_a0(double s)
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
      complex<double> temp = kernel(s, x[i]) * previous->interp_above(x[i]);
      sum += w[i] * temp;
    }
    return sum;
};

complex<double> inhomogeneity::integ_a0_a(double s)
{
  if (s < a0|| s > kinematics.pseudo_threshold())
  {
    cout << "integ_a0_a: Integration out of range! Quitting... \n";
    exit(1);
  }

  double wM[N_integ + 1], xM[N_integ + 1];
  gauleg(t_minus(s), sthPi - EPS, wM, xM, N_integ);

  double wP[N_integ + 1], xP[N_integ + 1];
  gauleg(sthPi - EPS, t_plus(s), wP, xP, N_integ);

  complex<double> sumP = 0., sumM = 0.;
  for(int i = 1; i < N_integ + 1; i++)
  {
    complex<double> tempM = kernel(s, xM[i]) * previous->interp_below(xM[i]);
    complex<double> tempP = kernel(s, xP[i]) * previous->interp_above(xP[i]);

    sumM += wM[i] * tempM;
    sumP += wP[i] * tempP;
  }

  return sumM + sumP;
};

complex<double> inhomogeneity::integ_a_b(double s)
{
  if (s < kinematics.pseudo_threshold() || s > kinematics.threshold())
  {
    cout << "integ_a_b: Integration out of range! Quitting... \n";
    exit(1);
  }

  previous->omega.set_ieps(0);
  double wM[N_integ + 1], xM[N_integ + 1];
  gauleg(t_minus(s), 2. * t_minus(kinematics.threshold()), wM, xM, N_integ);

  double wP[N_integ + 1], xP[N_integ + 1];
  gauleg(2. * t_plus(kinematics.threshold()), t_plus(s), wP, xP, N_integ);

  complex<double> sumP = 0., sumM = 0.;
  for(int i = 1; i < N_integ + 1; i++)
  {
    complex<double> tempM = kernel(s, xM[i]) * previous->omega.eval(xM[i]);
    complex<double> tempP = kernel(s, xP[i]) * previous->omega.eval(xP[i]);

    sumM += wM[i] * tempM;
    sumP += wP[i] * tempP;
  }

  return sumM + sumP;
};

complex<double> inhomogeneity::integ_b(double s)
{
  if (s < kinematics.threshold())
  {
    cout << "integ_s0_a0: Integration out of range! Quitting... \n";
    exit(1);
  }

  previous->omega.set_ieps(0);
  double w[N_integ + 1], x[N_integ + 1];
  gauleg(t_minus(s), t_plus(s), w, x, N_integ);

  complex<double> sum = 0.;
  for (int i = 1; i < N_integ + 1; i++)
  {
    complex<double> temp = kernel(s, x[i]) * previous->omega.eval(x[i]);
    sum += w[i] * temp;
  }
  return sum;
};

// ---------------------------------------------------------------------------
// Angular integral, F_hat
complex<double> inhomogeneity::operator()(double s)
{
  if (abs(s - sthPi) < EPS)
  {
    return 2. * previous->interp_above(t_minus(sthPi)) * pow((kinematics.pseudo_threshold() - sthPi), 1.5);
  }
  else if (s > sthPi && s < a0)
  {
    return integ_s0_a0(s);
  }
  else if (s >= a0 && s <= kinematics.pseudo_threshold())
  {
    return integ_a0_a(s);
  }
  else if (s > kinematics.pseudo_threshold() && s < kinematics.threshold())
  {
    return integ_a_b(s);
  }
  else if (s >= kinematics.threshold())
  {
    return integ_b(s);
  }
};
