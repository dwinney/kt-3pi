// Routines to fit a lineshape to a polynomial (z, theta) expansion around center of dalitz plot.
//
// Dependencies: dalitz.cpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "poly_param_fit.hpp"

// Polynomial expansion around the center of the dalitz plot.
// The bounds of integrations are easiest seen in terms of s and t (c.f. PDG Kinematics)
template <class T>
double poly_param_fit<T>::F_poly(double s, double t)
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
template <class T>
void poly_param_fit<T>::set_integration_points(int m)
{
  n = m;
  cout << "Number of integration points changed to " << n << "... \n";
  S_WG_GENERATED = false;
};

// Set the number of parameters in the polynomial expansion.
// Default is 2, i.e. alpha and beta up to order 3/2 in z.
template <class T>
void poly_param_fit<T>::set_params(int n, const double *par)
{
  n_params = n;
  switch (N_params())
   {
    case 4: delta = par[3];
    case 3: gamma = par[2];
    case 2: beta = par[1];
    case 1: alpha = par[0];
  };
};

// Polynomial expansion around the center of dalitz plot
template <class T>
void poly_param_fit<T>::print_params(int a)
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

// ---------------------------------------------------------------------------
// The gauleg function in aux_math.hpp indexes from 1 to N so thats why
// these wrapper functions exist. Also to put them into vectors such that the number of points can change.
template <class T>
void poly_param_fit<T>::generate_s_weights()
{
  double weights[N_int() + 1], abscissas[N_int() + 1];

  double smn = dalitz<T>::amp.smin() + dalitz<T>::offset;
  double smx = dalitz<T>::amp.smax() - dalitz<T>::offset;

  gauleg(smn , smx, abscissas, weights, N_int() + 1);

  s_wgt.clear(); s_abs.clear();

  for (int i = 1; i < N_int() + 1; i++)
  {
          s_wgt.push_back(weights[i]);
          s_abs.push_back(abscissas[i]);
  }
  cout << "Gaussian weights for s in the Dalitz Region generated with "
  << N_int() << " points... \n";
  S_WG_GENERATED = true;
};

// Same but now the bounds in the t variable depend on s.
template <class T>
void poly_param_fit<T>::generate_t_weights(vector<double> s)
{
  // Clear preexisting weights and abcsissas
  t_wgt.clear(); t_abs.clear();
  for (int i = 0; i < N_int(); i++)
  {
    double weights[N_int() + 1], abscissas[N_int() + 1];

    // add 0.001 to push values off the boundary where things may be singular
    double tmn = dalitz<T>::amp.tmin(s[i]) + dalitz<T>::offset;
    double tmx = dalitz<T>::amp.tmax(s[i]) - dalitz<T>::offset;

    gauleg(tmn, tmx, abscissas, weights, N_int() + 1);

    vector<double> t_wgt_temp, t_abs_temp;
    for (int j = 1; j < N_int() + 1; j++)
    {
        t_wgt_temp.push_back(weights[j]);
        t_abs_temp.push_back(abscissas[j]);
    }

    t_wgt.push_back(t_wgt_temp);
    t_abs.push_back(t_abs_temp);
  }
  cout << "Gaussian weights for t in the Dalitz Region generated with "
  << N_int() << " points... \n";
  T_WG_GENERATED = true;
};

template <class T>
void poly_param_fit<T>::generate_weights(){
  if (S_WG_GENERATED == false)
  {
    generate_s_weights();
    generate_t_weights(s_abs);
  }
  else if (T_WG_GENERATED == false)
  {
    generate_t_weights(s_abs);
  }
};

//-----------------------------------------------------------------------------
// The kinematic kernel that goes into the integral over the Dalitz region.
// For the omega case it is the kibble function, but this may be different in general.
template <class T>
double poly_param_fit<T>::kin_kernel(double s, double t)
{
  double temp1, temp2;
  temp1 = dalitz<T>::amp.Kibble(s,t);
  temp2 = dalitz<T>::amp.Kibble( dalitz<T>::amp.s_c(), dalitz<T>::amp.t_c()); // Normalized to the center of the Dalitz region
  return temp1 / temp2;
};

// Calculate the area of the physical Dalitz region
template <class T>
double poly_param_fit<T>::dalitz_area()
{
double t_sum, s_sum, s_i, t_ij;
generate_weights();

s_sum = 0.;
for (int i = 0; i < N_int(); i++)
{
  s_i = s_abs[i];
  t_sum = 0.;
  for (int j = 0; j < N_int(); j++)
  {
    t_ij = t_abs[i][j];
    t_sum +=  t_wgt[i][j] * t_ij;
  }
  s_sum += s_wgt[i] * t_sum;
};
return s_sum;
};

//-----------------------------------------------------------------------------
// Calculate chi_squared
// Normalized by the area of the Dalitz plot.
template <class T>
double poly_param_fit<T>::chi_squared(const double *par)
{
  // Update parameters at the fitting step:
  set_params(N_params(), par);

  complex<double> amp_ij;
  double t_sum, s_sum, s_i, t_ij, d_area;
  double ampsqr_ij, poly_ij, kern_ij;
  double tmp1, tmp2, tmp3;
  double chi2;

  generate_weights();

  d_area = dalitz_area();

  // Integrate over s
  s_sum = 0.;
  for (int i = 0; i < N_int(); i++)
  {
    s_i = s_abs[i];
    // Integrate over t
    t_sum = 0.;
    for (int j = 0; j < N_int(); j++)
    {
        t_ij = t_abs[i][j];
        amp_ij = dalitz<T>::amp(s_i, t_ij);
        ampsqr_ij = real(amp_ij * conj(amp_ij));
        poly_ij = F_poly(s_i, t_ij);
        kern_ij = kin_kernel(s_i, t_ij);

        tmp1 = ampsqr_ij - poly_ij;
        tmp2 = tmp1 * (kern_ij * kern_ij);
        tmp3 = tmp2 * tmp2;

        t_sum += t_wgt[i][j] * tmp3;
    };
    s_sum += s_wgt[i] * t_sum;
  };

  chi2 = s_sum / (pow(Norm, 4.) * d_area);
  return chi2;
};
// ---------------------------------------------------------------------------
// Minimize the above chi_squared by calling Minuit2 in the ROOT::Math::Minimizer class
template <class T>
void poly_param_fit<T>::fit_params()
{
  ROOT::Math::Minimizer* minuit = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
  minuit->SetMaxFunctionCalls(100000);
  minuit->SetTolerance(0.0001);
  minuit->SetPrintLevel(0);

  ROOT::Math::Functor fcn(this, &poly_param_fit::chi_squared, N_params());
  minuit->SetFunction(fcn);
  cout << "Minimizing with " << N_params() << " free parameters... \n";
  switch (N_params())
  {
    case 4: minuit->SetVariable(3, "delta", 0., .01);
    case 3: minuit->SetVariable(2,"gamma", 0., 1.);
    case 2: minuit->SetLowerLimitedVariable(1,"beta", 0., .01, 0.001);
    case 1: {minuit->SetLowerLimitedVariable(0,"alpha", 0., .01, 0.001); break;}
    default:
    {
      cout << "fit_params(): Invalid number of free parameters in Minuit2. Quitting... \n";
      exit(1);
    }
  }

  minuit->Minimize();
  cout << "Minimizing done. " << endl;
  cout << endl;

  double chi2 = minuit ->MinValue();
  cout << "sqrt(chi2) = " << sqrt(chi2) << endl;
  cout << endl;
  minuit->PrintResults();
};
