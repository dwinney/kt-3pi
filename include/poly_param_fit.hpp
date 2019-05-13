// Routines to fit a lineshape to a polynomial (z, theta) expansion around center of dalitz plot.
//
// Dependencies: dalitz.cpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _POLY_FIT_
#define _POLY_FIT_

#include <vector>
#include "dalitz.hpp"
#include "aux_math.hpp"
#include <iomanip>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using std::vector;
using std::setw;

template <class T>
class poly_param_fit : public dalitz<T>
{
protected:
  // Polynomial expansion parameters up to order 5/2 in z
  int n_params = 2;
  int N_params(){return n_params;};
  double Norm = 1.451, alpha = 0., beta = 0., gamma = 0., delta = 0.;
  double scale = 1.;

  // For integration
  int n = 60; // Number of integration points
  int N_int(){
    return n;
  };

  bool S_WG_GENERATED, T_WG_GENERATED;
  vector<double> s_wgt, s_abs;
  vector< vector<double>> t_wgt, t_abs;
  void generate_s_weights();
  void generate_t_weights(vector<double> s);
  void generate_weights();

// ---------------------------------------------------------------------------
public:
  poly_param_fit(T my_amp) : dalitz<T>(my_amp){};

  void set_params(int n, const double *par);
  void print_params(int a = 0);
  void set_integration_points(int n);

  double F_poly(double s, double t);

  double kin_kernel(double s, double t); // Kinematic Kernal in dalitz region integral
  double chi_squared(const double *par); // Chi-squared between input line-shape and polynomial.
  double dalitz_area();

  void fit_params();
// ---------------------------------------------------------------------------
};
// ---------------------------------------------------------------------------
#endif
