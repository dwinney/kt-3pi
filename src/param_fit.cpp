// Routines to fit a lineshape to a polynomial (z, theta) expansion around center of dalitz plot.
//
// Dependencies: dalitz.cpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "param_fit.hpp"

// Polynomial expansion around the center of the dalitz plot.
// The bounds of integrations are easiest seen in terms of s and t (c.f. PDG Kinematics)
template <class T>
double param_fit<T>::F_poly(double s, double t)
{
  double zs, thetas, poly;
  zs = dalitz<T>::amp.z(s,t);
  thetas = dalitz<T>::amp.theta(s,t);

  poly = 1.0
        + 2. * alpha * scale * zs
        + 2. * beta * scale * std::pow(zs, 1.5) * std::sin(3.*thetas)
        + 2. * gamma * scale *  zs*zs
        + 2. * delta * scale * std::pow(zs, 2.5) * std::sin(3.*thetas);
  return Norm * Norm * poly;
};

// ---------------------------------------------------------------------------
// General set and print functions

// Set the number of parameters in the polynomial expansion.
// Default is 2, i.e. alpha and beta up to order 3/2 in z.
template <class T>
void param_fit<T>::set_params(int n, const double *par)
{
  switch (n)
   {
    case 4: delta = par[4];
    case 3: gamma = par[3];
    case 2: beta = par[2];
    case 1: alpha = par[1];
    case 0: Norm = par[0]; break;
  };
};

// Polynomial expansion around the center of dalitz plot
template <class T>
void param_fit<T>::print_params(int a)
{
  switch (a)
  {
  case 0: {
      cout << "Printing Dalitz Plot parameters... \n";
      cout << "---------------------------------------- \n";
      cout << std::left <<  setw(20) << "normalization:" << setw(20) << Norm << "\n";
      cout << std::left <<  setw(20) << "alpha:" << setw(20) << alpha << "\n";
      cout << std::left <<  setw(20) << "beta:" << setw(20) << beta << "\n";
      cout << std::left <<  setw(20) << "gamma:" << setw(20) <<  gamma << "\n";
      cout << std::left <<  setw(20) << "delta:" << setw(20) << delta << "\n";
      cout << "---------------------------------------- \n";
      cout << "\n";
    };
  };
};

//-----------------------------------------------------------------------------
// The kinematic kernel that goes into the integral over the Dalitz region.
// For the omega case it is the kibble function, but this may be different in general.
template <class T>
double param_fit<T>::kin_kernel(double s, double t)
{
  complex<double> temp1, temp2;
  temp1 = dalitz<T>::amp.Kibble(s,t);
  temp2 = dalitz<T>::amp.Kibble( dalitz<T>::amp.s_c(), dalitz<T>::amp.t_c()); // Normalized to the center of the Dalitz region
  return abs(temp1 / temp2);
};


//-----------------------------------------------------------------------------
// Calculate chi_squared
// Normalized by the area of the Dalitz plot.
template <class T>
double param_fit<T>::chi_squared(const double *par)
{
  // Update parameters at the fitting step:
  set_params(n_params, par);

  double t_sum, s_sum, s_i, d_area;
  double chi2;

  dalitz<T>::generate_weights();

  d_area = dalitz<T>::dalitz_area();

  // Integrate over s
  s_sum = 0.;
  for (int i = 0; i < dalitz<T>::N_int(); i++)
  {
    s_i = dalitz<T>::s_abs[i];
    // Integrate over t
    t_sum = 0.;
    for (int j = 0; j < dalitz<T>::N_int(); j++)
    {
        complex<double> amp_ij;
        double t_ij, ampsqr_ij, poly_ij, kern_ij;

        t_ij = dalitz<T>::t_abs[i][j];
        amp_ij = dalitz<T>::amp(s_i, t_ij);
        ampsqr_ij = abs(amp_ij * amp_ij);
        poly_ij = F_poly(s_i, t_ij);
        kern_ij = kin_kernel(s_i, t_ij);

        double tmp1, tmp2, tmp3;
        tmp1 = (ampsqr_ij - poly_ij);
        tmp2 = tmp1 * kern_ij;
        tmp3 = tmp2 * tmp2;

        t_sum += dalitz<T>::t_wgt[i][j] * tmp3;
    };
    s_sum += dalitz<T>::s_wgt[i] * t_sum;
  };

  chi2 = s_sum / d_area;

  return chi2;
};
// ---------------------------------------------------------------------------
// Minimize the above chi_squared by calling Minuit2 in the ROOT::Math::Minimizer class
template <class T>
void param_fit<T>::extract_params(double N)
{
  ROOT::Math::Minimizer* minuit = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
  minuit->SetMaxFunctionCalls(10000000);
  minuit->SetTolerance(0.001);
  minuit->SetPrintLevel(0);

  n_params = N;
  ROOT::Math::Functor fcn(this, &param_fit::chi_squared, N);
  minuit->SetFunction(fcn);
  cout << "param_fit: Minimizing with " << N << " free parameters... \n";

  minuit->SetVariable(0,"normalization", 1., .01);
  if (N >= 1)
  {
    minuit->SetVariable(1,"alpha", 0., .01);
    if (N >= 2)
    {
      minuit->SetVariable(2,"beta", 0., .01);
      if (N >= 3)
      {
        minuit->SetVariable(3,"gamma", 0., 1.);
        if (N == 4)
        {
          minuit->SetVariable(4, "delta", 0., .01);
          if (N > 4)
          {
            cout << "fit_params(): Invalid number of free parameters in Minuit2. Quitting... \n";
            exit(1);
          }
        }
      }
    }
  }
else
  {
    cout << "fit_params(): Invalid number of free parameters in Minuit2. Quitting... \n";
    exit(1);
  }

  minuit->Minimize();
  cout << "Minimizing done. " << endl;
  cout << endl;

  double chi2 = minuit ->MinValue();
  cout << "sqrt(chi2) = " << sqrt(chi2) << endl;
  cout << endl;

  print_params();

};
