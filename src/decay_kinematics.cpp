// Generic kinematics class for the reaction of J^PC -> 3 pi.
// The reaction is defined by the mass and the quantum numbers of the decaying particle.
// This allows passing all relevant kinematic quantities easily into other classes.
//
// Dependencies: constants.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "decay_kinematics.hpp"

// Mandelstam t in terms of s-channel CM frame variables
complex<double> decay_kinematics::t_man(complex<double> s, complex<double> zs)
{
  complex<double> ks = sqrt(Kallen_pi(s) * Kallen_x(s)) / (4. * s);
  return (3. * mPi * mPi + mDec * mDec  - s) / 2. + 2. * ks * zs;
};

// Mandelstam u in terms of the other two invariant variables
complex<double> decay_kinematics::u_man(complex<double> s, complex<double> t)
{
 return 3.*mPi*mPi + mDec*mDec - s - t;
};

// Lorentz Covariant Kibble function
complex<double> decay_kinematics::Kibble(complex<double> s, complex<double> t)
{
    complex<double> u = u_man(s,t);
    complex<double> lambda = (s * t * u) - mPi*mPi * (mDec*mDec - mPi*mPi)*(mDec*mDec - mPi*mPi);
    return xr * lambda;
};

complex<double> decay_kinematics::Kallen(complex<double> x, complex<double> y, complex<double> z)
{
  return x*x + y*y + z*z - 2.*(x*y + x*z + z*y);
};

//-----------------------------------------------------------------------------
// Scattering angles in the s, t, and u center-of-mass frames
double decay_kinematics::z_s(double s, double t)
{
  complex<double> ks = sqrt(Kallen_x(s) * Kallen_pi(s)) / s;
  complex<double> result =  (t - u_man(s,t)) / ks;

  if (abs(imag(result)) > 0.0001)
  {
    cout << "z_s: error arguments are complex. Quitting..." << endl;
    exit(1);
  }

  return real(result);
};

double decay_kinematics::z_t(double s, double t)
{
  return z_s(t , s);
};

double decay_kinematics::z_u(double s, double t)
{
  double u = real(u_man(s,t));
  return z_s(u, t);
};

// ---------------------------------------------------------------------------
// Momenta and Energies in the Center of Momentum Frame
complex<double> decay_kinematics::com_E2(double s)
{
  return sqrt(s * xr) / 2.;
};

complex<double> decay_kinematics::com_E3(double s)
{
  return (mDec*mDec - mPi*mPi - s * xr) / (2.*sqrt(s * xr));
};

complex<double> decay_kinematics::com_P2(double s)
{
  complex<double> temp = com_E2(s)*com_E2(s) - mPi*mPi;
  return sqrt(temp);
};

complex<double> decay_kinematics::com_P3(double s)
{
  complex<double> temp = com_E3(s)*com_E3(s) - mPi*mPi;
  return sqrt(temp);
};

double decay_kinematics::tmin(double s)
{
  return real(pow(com_E2(s) + com_E3(s), 2.) - pow(com_P2(s) + com_P3(s), 2.));
};

double decay_kinematics::tmax(double s)
{
  return real(pow(com_E2(s) + com_E3(s), 2.) - pow(com_P2(s) - com_P3(s), 2.));
};

//-----------------------------------------------------------------------------
// Angular functions

// Outputs the Wigner little-d function appropriate for the decay into three scalar,
// d^j_{lambda, 0}(z)
double decay_kinematics::d_func0(int j, int l, double z)
{
  double prefactor = TMath::Factorial(j - l) / TMath::Factorial(j + l);

  double legendre;
  if (l >= 0)
  {
    legendre = ROOT::Math::assoc_legendre(j, l, z);
  }
  else
  {
    legendre = ROOT::Math::assoc_legendre(j, -l, z);
    legendre *= pow(-1., double(l));
    legendre *= TMath::Factorial(j - l) / TMath::Factorial(j + l);
  }

  return sqrt(2.) * sqrt(prefactor) * legendre;
};

// Outputs d_hat the kinematic-singularity-free d_function
double decay_kinematics::d_hat(int j, int l, double z)
{
  double d_lam = d_func0(j, l, z);
  double half_angle_factor = pow((1. - z*z), abs(double(l))/2.);

  return d_lam / half_angle_factor;
};

complex<double> decay_kinematics::d_hat(int j, int l, complex<double> z)
{
  if (l == 1)
  {
    switch (j)
    {
    case 1: return 1.;
    case 3: return 3. * z;
    default: cout << "d_hat: complex legendres not coded that far yet... \n";
              exit(0);
    }
  }
  else
  {
    cout << "d_hat: complex legendres not coded that far yet... \n";
    exit(0);
  }
};

//-----------------------------------------------------------------------------
// Angular momentum barrier factor
complex<double> decay_kinematics::barrier_factor(int j, int lam, complex<double> s)
{
  double dlam = double(lam), dj = double(j);

  // Angular momentum barrier factor
  complex<double> temp2 = sqrt(Kallen_x(s) * Kallen_pi(s));
  return pow(temp2, dj - abs(dlam));
};

//-----------------------------------------------------------------------------
// Kinematic Singularities of Helicity amplitudes
// helicity projection : Lambda
// spin projection : j > 0
complex<double> decay_kinematics::K_jlam(int j, int lam, complex<double> s, complex<double> zs)
{
  if (j < std::abs(lam))
  {
    cout << "K_jlam: j < |lambda| is not allowed. Quitting..." << endl;
    exit(1);
  }

  double dlam = double(lam), dj = double(j);

  // Half angle factors and K factors
  complex<double> temp = Kibble(s, t_man(s, zs)) / 4.;
  complex<double> result = pow(sqrt(temp), abs(dlam));

  // LS - lambda little group factor
  int Yx = qn_J - (1 + get_naturality())/2;
  complex<double> temp3 = sqrt(Kallen_x(s));
  result *= pow(temp3, abs(dj - double(Yx)) - dj);

  return result;
};

// ---------------------------------------------------------------------------
// Lorentz Invariant dimensionless parameters in terms of Mandelstam variables
double decay_kinematics::x(double s, double t)
{
  double temp;
  temp = sqrt(3.) * real((t - u_man(s,t)));
  return temp * d_norm();
};

double decay_kinematics::y(double s, double t)
{
  double temp;
  temp = 3. * (s_c() - s);
  return temp * d_norm();
};

//-----------------------------------------------------------------------------
// Lorentz Invariant dimensionless parameters in terms of polar variables
double decay_kinematics::x_polar(double r, double phi)
{
  return sqrt(r) * cos(phi);
};

double decay_kinematics::y_polar(double r, double phi)
{
  return sqrt(r) * sin(phi);
};
//-----------------------------------------------------------------------------
// Inverted to get r and phi in terms of s and t
double decay_kinematics::r(double s, double t)
{
    double X, Y;
    X = x(s,t);
    Y = y(s,t);
    return X * X + Y * Y;
};

double decay_kinematics::phi(double s, double t)
{
  double result = atan2(y(s,t), x(s,t));
  if (result < 0.0)
  {
    result += 2. * M_PI; // Map [-pi,pi] to [0, 2pi]
  }
  return result;
};
