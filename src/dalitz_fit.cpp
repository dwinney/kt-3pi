// Routines to fit a lineshape to another lineshape by minimizing a chi2-type integral over
// the entire dalitz region
//
// Dependencies: dalitz.cpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dalitz_fit.hpp"

//-----------------------------------------------------------------------------
// The kinematic kernel that goes into the integral over the Dalitz region.
// For the omega case it is the kibble function, but this may be different in general.
template <class T, class F>
double dalitz_fit<T, F>::kin_kernel(double s, double t)
{
  complex<double> temp1, temp2;
  temp1 = dalitz<T>::amp->kinematics.Kibble(s,t);
  temp2 = dalitz<T>::amp->kinematics.Kibble( dalitz<T>::amp->kinematics.s_c(), dalitz<T>::amp->kinematics.t_c()); // Normalized to the center of the Dalitz region
  return abs(temp1 / temp2);
};


//-----------------------------------------------------------------------------
// Calculate chi_squared
// Normalized by the area of the Dalitz plot.
template <class T, class F>
double dalitz_fit<T, F>::chi_squared(const double *par)
{
  dalitz<T>::generate_weights();

  // Update parameters in polynomial expansion at the fitting step
  fit_amp->set_params(n_params, par);

  // Integrate over s
  double s_sum = 0.;
  for (int i = 0; i < dalitz<T>::N_int(); i++)
  {
    double s_i = dalitz<T>::s_abs[i];
    // Integrate over t
    double t_sum = 0.;
    for (int j = 0; j < dalitz<T>::N_int(); j++)
    {
        complex<double> amp_ij, fit_ij;
        double t_ij, ampsqr_ij, fitsqr_ij, kern_ij;

        t_ij = dalitz<T>::t_abs[i][j];

        amp_ij = dalitz<T>::amp->eval(s_i, t_ij);
        ampsqr_ij = abs(amp_ij * amp_ij);

        fit_ij = fit_amp->eval(s_i, t_ij);
        fitsqr_ij = abs(fit_ij * fit_ij);

        kern_ij = kin_kernel(s_i, t_ij);

        double tmp, tmp3;
        tmp = (ampsqr_ij - fitsqr_ij);
        tmp *= kern_ij * kern_ij;
        tmp /= dalitz<T>::amp->error_func(s_i, t_ij);
        tmp /= par[0] * par[0];

        tmp3 = tmp * tmp;

        t_sum += dalitz<T>::t_wgt[i][j] * tmp3;
    };
    s_sum += dalitz<T>::s_wgt[i] * t_sum;
  };

  double chi2 = s_sum / dalitz<T>::dalitz_area();;

  return chi2;
};

// ---------------------------------------------------------------------------
// Minimize the above chi_squared by calling Minuit2 in the ROOT::Math::Minimizer class
template <class T, class F>
void dalitz_fit<T, F>::extract_params(int eN)
{
  ROOT::Math::Minimizer* minuit = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
  minuit->SetMaxFunctionCalls(10000000);
  minuit->SetTolerance(0.001);
  minuit->SetPrintLevel(nError);

  n_params = eN;
  ROOT::Math::Functor fcn(this, &dalitz_fit::chi_squared, n_params);
  minuit->SetFunction(fcn);
  cout << "dalitz_fit: Fitting";

  if (dalitz<T>::amp->kinematics.get_ampName() != "")
  {
    cout << " " + dalitz<T>::amp->kinematics.get_ampName();
  }
   cout << " with " << n_params << " free parameters... \n";

  for (int a = 0; a < n_params ; a++)
  {
    string var_name = "par[" + std::to_string(a) + "]";
    minuit->SetVariable(a, var_name.c_str(), 1., 0.1);
  }

  minuit->Minimize();
  cout << "Minimizing done. " << endl;
  cout << endl;

  double chi2 = minuit ->MinValue();
  cout << "sqrt(chi2) = " << sqrt(chi2) * 1.e3 << endl;
  cout << endl;

  fit_amp->set_params(n_params, minuit->X());

};
