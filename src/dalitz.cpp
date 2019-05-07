// Classes for Dalitz plot generation and fitting to polynomial expansion with dalitz plot parameters.
//
// Dependencies: amp.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dalitz.hpp"

//-----------------------------------------------------------------------------
// Lorentz Invariant dimensionless parameters in terms of Mandelstam variables
double dalitz::x(double s, double t)
{
  temp1 = sqrt(3.) * (t - u(s,t) );
  temp2 = 2.*mDec* (mDec - 3.*mPi);

  return temp1 / temp2;
};

double dalitz::y(double s, double t)
{
  temp1 = 3. * (s_c - s);
  temp2 = 2.* mDec * (mDec - 3.*mPi);

  return temp1 / temp2;
};

// Lorentz Invariant dimensionless parameters in terms of polar variables
double dalitz::x_polar(double z, double theta)
{
  return sqrt(z) * cos(theta);
};

double dalitz::y_polar(double z, double theta)
{
  return sqrt(z) * sin(theta);
};

// Inverted to get z and theta in terms of s and t
double dalitz::z(double s, double t)
{
    temp1 = x(s,t);
    temp2 = y(s,t);

    return temp1 * temp1 + temp2 * temp2;
};

double dalitz::theta(double s, double t)
{
  temp1 = x(s,t);
  temp2 = y(s,t);

  return atan2(temp2, temp1);
};
//-----------------------------------------------------------------------------

// Doubly Differential Decay Width
double dalitz::d2Gamma(double s, double t)
{
  temp1 = 96. * (mDec*mDec*mDec) * pow(2.* M_PI, 3);
  amp_ij = abs(amp(s,t));
  return amp_ij * amp_ij / temp1;
};

//-----------------------------------------------------------------------------
// Extraction of Dalitz Plot parameters by fitting
//-----------------------------------------------------------------------------
// Polynomial expansion around the center of the dalitz plot.
// The bounds of integrations are easiest seen in terms of s and t (c.f. PDG Kinematics)
double dalitz::F_poly(double s, double t)
{
  temp1 = z(s,t);
  temp2 = theta(s,t);

  temp3 = 1.0
        + 2. * alpha * scale * temp1
        + 2. * beta * scale * pow(temp1, 1.5) * sin(3.*temp2)
        + 2. * gamma * scale *  temp1*temp1
        + 2. * delta * scale * pow(temp1, 2.5) * sin(3.*temp2);
  return Norm * Norm * temp3;
};
//-----------------------------------------------------------------------------
// The gauleg function in aux_math.hpp indexes from 1 to N so thats why
// these wrapper functions exist. Also to put them into vectors such that the number of points can change.
void dalitz::generate_s_weights()
{
  double weights[N_int + 1], abscissas[N_int + 1];
  gauleg(smin, smax, abscissas, weights, N_int + 1);

  s_wgt.clear(); s_abs.clear();

  for (int i = 1; i < N_int + 1; i++)
  {
          s_wgt.push_back(weights[i]);
          s_abs.push_back(abscissas[i]);
  }

  S_WG_GENERATED = true;
};

void dalitz::generate_t_weights(double s)
{
  double weights[N_int + 1], abscissas[N_int + 1];
  gauleg(tmin(s), tmax(s), abscissas, weights, N_int + 1);

  t_wgt.clear(); t_abs.clear();

  for (int i = 1; i < N_int + 1; i++)
  {
          t_wgt.push_back(weights[i]);
          t_abs.push_back(abscissas[i]);
  }
};
//-----------------------------------------------------------------------------

// The kinematic kernal that goes into the integral over the Dalitz region.
// For the omega case it is the kibble function, but this may be different in general.
double dalitz::kin_kernal(double s, double t)
{
  temp1 = Kibble(s,t);
  temp2 = Kibble(0., 0.);

  return temp1 / temp2;
};

// Calculate chi_squared
double dalitz::chi_squared(double s, double t)
{
  if (S_WG_GENERATED == false)
  {
    generate_s_weights();
  };

// Integrate over s
for (int i = 0; i < N_int; i++)
{
  s_i = s_abs[i];
  generate_t_weights(s_i);

  // Integrate over t
  t_sum = 0;
  for (int j = 0; j < N_int; j++)
  {
      t_j = t_abs[j];
      amp_ij = abs(amp(s_i, t_j));
      poly_ij = F_poly(s_i, t_j);

      temp1 = amp_ij * amp_ij - poly_ij * poly_ij;
      temp2 = kin_kernal(s_i, t_j) * temp1 * temp1;

      t_sum += t_wgt[j] * temp2;
      d_area += t_wgt[j] * t_j; // The area of the dalitz plot to normalize the chi2 value.
  };
  s_sum += s_wgt[i] * t_sum;
  d_area += s_wgt[i] * s_i;
};
temp3 = s_sum / (Norm * Norm * d_area);
return temp3;
};
