// These are methods for evaluating the KT angular (t) integral (i.e. the angular projection of cross subchannel
// isobars). This is seperated to allow the possibility of testing different methods of evaluating dispersion (s) integral.
//
// Dependencies: kt_iteration.hpp, decay_kinematics.hpp, aux_math.hpp, omnes.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "kt_ang_integral.hpp"

// ---------------------------------------------------------------------------
// Kacser function k(s)
complex<double> angular_integral::k(complex<double> s)
{
  complex<double> temp1 = sqrt(xr * (s - sthPi + EPS) / s);
  temp1 *= sqrt(xr * (b - s)) * sqrt(xr * (a - s));

  return temp1;
};

// ---------------------------------------------------------------------------
// COM Scattering angle in the s-channel scattering subchannel
// Complex in general because of the path of integration needed by analytic continuation
complex<double> angular_integral::sine_sqr(double s, complex<double> t)
{
  complex<double> zs = 2. * t + s - mDec * mDec - 3. * mPi * mPi;
  zs /= k(s);

  return 1. - zs * zs;
};

// ---------------------------------------------------------------------------
// Bounds of integration
complex<double> angular_integral::t_minus(double s)
{
    return (mDec * mDec + 3. * mPi * mPi - s) / 2. - k(s) / 2.;
};

complex<double> angular_integral::t_plus(double s)
{
    return (mDec * mDec + 3. * mPi * mPi - s) / 2. + k(s) / 2.;
};

// ---------------------------------------------------------------------------
// Kinematic kernel inside the angular integral. In general this should depend
// on the helicity and spin of the isobar.
//
// Right now it only has the omega case of Lambda = 1, j = 1, I = I
complex<double> angular_integral::kernel(int j, int lam, int jp, int lamp, double s, complex<double> t)
{
  complex<double> zs = 2. * t + s - mDec * mDec - 3. * mPi * mPi;
  zs /= k(s);

  complex<double> zt = 2. * s + t - mDec * mDec - 3. * mPi * mPi;
  zt /= k(t);

  complex<double> result = (2. * jp + 1.);
  result *= pow(sine_sqr(s, t), double(lam));

  // Ratio of kibble functions is easy
  result *= pow(kinematics.Kibble(s,t) / 4., double(lam - lamp)/2.);

  // Ratio of Little Group factors
  int Yx = kinematics.get_totalSpin() - (1 + kinematics.get_naturality())/2;
  int Yxj = abs(j - Yx) - j;
  int Yxjp = abs(jp - Yx) - jp;
  result *= pow(kinematics.Kallen_x(t), double(Yxjp)/2.) / pow(kinematics.Kallen_x(s), double(Yxj)/2.);

  // Angular function in t- channel
  result *= kinematics.barrier_factor(jp, lamp, t);
  result *= kinematics.d_hat(jp, lamp, zt);

  // extra factor of k(s) from the change of variables from z_s to t (i.e. Jacobian)
  return result / k(s);
};

// ---------------------------------------------------------------------------
// Integration path in the complex plane

// Both limits are real and above the unitarity cut
complex<double> angular_integral::integ_sthPi_a0(int j, int n, double s)
{
    if (s <= sthPi || s >= a0)
    {
      cout << "integ_s0_a0: Integration out of range! Quitting... \n";
      exit(1);
    }

    double w[N_integ + 1], x[N_integ + 1];
    gauleg(real(t_minus(s)), real(t_plus(s)), x, w, N_integ);

    complex<double> sum = 0.;

    for (int jp = 0; 2*jp+1 <= options.max_spin; jp++)
    {
      for (int i = 1; i < N_integ + 1; i++)
      {
        complex<double> temp = kernel(2*j+1, 1, 2*jp+1, 1, s, x[i]) * previous->isobars[jp].subtractions[n].interp_above(x[i]);
        sum += w[i] * temp;
      }
    }
    return sum;
};

// both limits are purely real but one is above the other below
complex<double> angular_integral::integ_a0_a(int j, int n, double s)
{
  if (s < a0|| s > a )
  {
    cout << "integ_a0_a: Integration out of range! Quitting... \n";
    exit(1);
  }

  complex<double> sumM, sumP;

  for (int jp = 0; 2*jp+1 <= options.max_spin; jp++)
  {
    if (std::abs(t_minus(s) - (sthPi + EPS)) > 0.001)
    {
      double wM[N_integ + 1], xM[N_integ + 1];
      gauleg(real(t_minus(s)), sthPi + EPS, xM, wM, N_integ);

      for(int i = 1; i < N_integ + 1; i++)
      {
        complex<double> tempM = kernel(2*j+1, 1, 2*jp+1, 1, s, xM[i]) * previous->isobars[jp].subtractions[n].interp_below(xM[i]);
        sumM += wM[i] * tempM;
      }
    }

    if (std::abs(t_plus(s) - (sthPi + EPS)) > 0.001)
    {
      double wP[N_integ + 1], xP[N_integ + 1];
      gauleg(sthPi + EPS, real(t_plus(s)), xP, wP, N_integ);

      for(int i = 1; i < N_integ + 1; i++)
      {
        complex<double> tempP = kernel(2*j+1, 1, 2*jp+1, 1, s, xP[i]) * previous->isobars[jp].subtractions[n].interp_above(xP[i]);
        sumP += wP[i] * tempP;
      }
    }
  }

  return sumM + sumP;
};

// this is the unphysical region, the bounds of integration are complex to avoid singularities
complex<double> angular_integral::integ_a_b(int j, int n, double s)
{
  if (s < a || s > b)
  {
    cout << "integ_a_b: Integration out of range! Quitting... \n";
    exit(1);
  }

  double w[N_integ + 1], x[N_integ + 1];
  gauleg(0., 1., x, w, N_integ);

  complex<double> sumM = 0.;
  complex<double> sumP = 0.;

  for (int jp = 0; 2*jp+1 <= options.max_spin; jp++)
  {
    for(int i = 1; i < N_integ + 1; i++)
    {
      complex<double> z1_i = (1. - x[i]) * t_minus(s) + x[i] * (2. * t_minus(b));
      complex<double> tempM = kernel(2*j+1, 1, 2*jp+1, 1, s, z1_i) * previous->isobars[jp].omega(z1_i, 0) * sub_poly(n, z1_i, 0);

      tempM *= 2. * t_minus(b) - t_minus(s); // Jacobian
      sumM +=  w[i] * tempM;
    }


    for(int i = 1; i < N_integ + 1; i++)
    {
      complex<double> z2_i = (1. - x[i]) *  (2. * t_plus(b)) + x[i] * t_plus(s);
      complex<double> tempP = kernel(2*j+1, 1, 2*jp+1, 1, s, z2_i) * previous->isobars[jp].omega(z2_i, 0) * sub_poly(n, z2_i, 0);

      tempP *= t_plus(s) - 2. * t_plus(b);
      sumP +=  w[i] * tempP;
    }
  }

  return (sumP + sumM);
};

// t-channel scattering region, limits are real again
complex<double> angular_integral::integ_b(int j, int n, double s)
{
  if (s < b)
  {
    cout << "integ_b: Integration out of range! Quitting... \n";
    exit(1);
  }

  double w[N_integ + 1], x[N_integ + 1];
  gauleg(real(t_minus(s)), real(t_plus(s)), x, w, N_integ);

  complex<double> sum = 0.;

  for (int jp = 0; 2*jp+1 <= options.max_spin; jp++)
  {
    for (int i = 1; i < N_integ + 1; i++)
    {
      complex<double> temp = kernel(2*j+1, 1, 2*jp+1, 1, s, x[i]) * previous->isobars[jp].omega(x[i], 0) * sub_poly(n, x[i], 0);
      sum += w[i] * temp;
    }
  }

  return sum;
};

// ---------------------------------------------------------------------------
// Angular integral, F_hat.
// This is the discontinuity or the inhomogeneity of the isobar
complex<double> angular_integral::operator () (int j, int n, double s)
{
  complex<double> integ;

  if (s > sthPi && s < a0) {
    integ = integ_sthPi_a0(j, n, s);
  }

  else if (s >= a0 && s <= a)  {
    integ = integ_a0_a(j, n, s);
  }

  else if (s > a && s < b)  {
    integ = integ_a_b(j, n, s);
  }

  else if (s >= b)  {
    integ = integ_b(j, n, s);
  }

  else  {
    cout << "AAAAA wtf hella error" << endl;
    exit(1);
  }

  return integ;
};

// ---------------------------------------------------------------------------
// Way to pass the current iteration to the angular_integral
void angular_integral::pass_iteration(iteration * prev)
{
    previous = NULL; //reset the pointer for posterity

    previous = prev;
};
