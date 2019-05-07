// Classes for Dalitz plot generation and fitting to polynomial expansion with dalitz plot parameters.
//
// Dependencies: amp.hpp, aux_math.hpp, ROOT
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _DALITZ_
#define _DALITZ_

#include "amp.hpp"
#include "aux_math.hpp"
#include <cmath>
#include <vector>
// #include <TFitter.h>
// #include <TMinuit.h>
// #include <TROOT.h>

using std::vector;

//-----------------------------------------------------------------------------
//  Define a Dalitz plot object with a pointer to the amplitude.
//  The amp should be a complex<double> and a function of 2 doubles, Mandelstam s and t.
//
//  complex<double> (*ptr_amp) (double, double);
//  ptr_func = my_amp;
//
//  dalitz my_dalitz(ptr_amp);
//-----------------------------------------------------------------------------

class dalitz : public amplitude
{
protected:
complex<double> (*amp) (double, double);
double s_c = (mDec * mDec + 3. * mPi*mPi) / 3.;
double scale = 1.0e-3;

//Temporary variables
double temp1, temp2, temp3, amp_ij, poly_ij;
double t_sum, s_sum, s_i, t_j, d_area;

// Polynomial expansion parameters
double Norm = 1., alpha, beta, gamma, delta;

// Gaussian-Legendre Weights and Absiccas in s
bool S_WG_GENERATED;
double chi2;
vector<double> s_wgt, s_abs, t_wgt, t_abs;
int N_int = 60; // Number of integration points
int N_params = 2;

//-----------------------------------------------------------------------------
public:
  dalitz()
  {
    cout << " Need pointer to amplitude function in dalitz object declaration. Quitting... \n";
    std::exit(1);
  };
  dalitz(complex<double> (*my_amp) (double, double))
  {
    amp = my_amp;
  };

//-----------------------------------------------------------------------------
  // Lorentz Invariant dimensionless parameters
  double x(double s, double t);
  double y(double s, double t);

  double x_polar(double z, double theta);
  double y_polar(double z, double theta);

  double z(double s, double t);
  double theta(double s, double t);

  // Double differential cross section
  double d2Gamma(double s, double t);

//-----------------------------------------------------------------------------
// Extraction of Dalitz Plot parameters by fitting
//-----------------------------------------------------------------------------
  void set_integration_points(int n)
  {
    N_int = n;
    S_WG_GENERATED = false;
  };

  // Set the number of parameters in the polynomial expansion.
  // Default is 2, i.e. alpha and beta up to order 3/2 in z.
  void set_N_params(int n)
  {
    N_params = n;
  };

  void generate_s_weights();
  void generate_t_weights(double s);

  // Polynomial expansion around the center of dalitz plot
  double F_poly(double s, double t);
  void print_params()
  {
    cout << " Printing Dalitz Plot parameters... \n";
    cout << " ---------------------------------------- \n";
    cout << " normalization :  \t \t" << Norm << "\n";
    cout << " alpha : \t \t" << alpha << "\n";
    cout << " beta : \t \t " << beta << "\n";
    cout << " gamma : \t \t " << gamma << "\n";
    cout << " delta : \t \t" << delta << "\n";
    cout << " ---------------------------------------- \n";
    cout << "\n";
  };

  double kin_kernal(double s, double t); // Kinematic Kernal in dalitz region integral
  double chi_squared(); // Chi-squared between input line-shape and polynomial.

  // void root_WRAPPER(int& npar, double* g, double& result, double *par, int flag);
  // void fit_params();
};

//-----------------------------------------------------------------------------
#endif
