// General routines for Dalitz plot generation
//
// Dependencies: constants.pp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dalitz.hpp"

//-----------------------------------------------------------------------------
// Lorentz Invariant dimensionless parameters in terms of Mandelstam variables
template <class T>
double dalitz<T>::x(double s, double t)
{
  double temp;
  temp = sqrt(3.) * (t - amp.u(s,t));
  return temp * d_norm();
};

template <class T>
double dalitz<T>::y(double s, double t)
{
  double temp;
  temp = 3. * (s_c() - s);
  return temp * d_norm();
};
//-----------------------------------------------------------------------------
// Lorentz Invariant dimensionless parameters in terms of polar variables
template <class T>
double dalitz<T>::x_polar(double z, double theta)
{
  return sqrt(z) * cos(theta);
};

template <class T>
double dalitz<T>::y_polar(double z, double theta)
{
  return sqrt(z) * sin(theta);
};
//-----------------------------------------------------------------------------
// Inverted to get z and theta in terms of s and t
template <class T>
double dalitz<T>::z(double s, double t)
{
    double tmp1, tmp2;
    tmp1 = x(s,t);
    tmp2 = y(s,t);

    return tmp1 * tmp1 + tmp2 * tmp2;
};

template <class T>
double dalitz<T>::theta(double s, double t)
{
  double temp1, temp2, temp3;
  temp1 = x(s,t);
  temp2 = y(s,t);

  temp3 = atan2(temp2, temp1);
  if (temp3 < 0.0)
  {
    temp3 += 2. * M_PI; // Map [-pi,pi] to [0, 2pi]
  }
  return temp3;
};
//-----------------------------------------------------------------------------

// Doubly Differential Decay Width
template <class T>
double dalitz<T>::d2Gamma(double s, double t)
{
  complex<double> f;
  double temp1;
  temp1 = 96. * pow(2.* M_PI * amp.mDec, 3);
  f = amp(s,t);
  // cout << std::left << setw(15) << s << setw(15) << t << setw(15) << abs(f*f) * Kibble(s,t) << endl;
  return abs( f * f) * amp.Kibble(s,t) / temp1;
};

//-----------------------------------------------------------------------------
// Quick function to print out a txt file with the Dalitz plot.
template <class T>
void dalitz<T>::plot(const char * n)
{
  std::ofstream output;
  output.open(n);

  double s_step = (amp.smax() - amp.smin() - offset )/100.;
  double t_step;

  double si, tij, gam;

  for (int i = 0; i < 100; i++)
  {
    si = amp.smin() + offset + double(i) * s_step;
    for(int j = 0; j < 100; j++)
    {
      t_step = (amp.tmax(si) - offset - amp.tmin(si))/100.;
      tij = amp.tmin(si) + offset  + double(j) * t_step;

      gam = d2Gamma(si, tij);

      output << std::left << setw(15) << x(si, tij) << setw(15) << y(si, tij) << setw << 10* gam << endl;
    }
  }
output.close();
};
