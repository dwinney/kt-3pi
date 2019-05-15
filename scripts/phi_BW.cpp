#include "poly_param_fit.cpp"
#include "dalitz.cpp"
#include "breit_wigner.hpp"
#include "test_poly.hpp"

#include <iostream>
#include <complex>


using std::complex;
using std::cout;
using std::endl;

int main()
{

  // Theoretical amplitude is a simple Breit-Wigner for the rho
  breit_wigner_KLOE rho(.770, .150);
  rho.set_Mass(1.02); // Phi decay

  // Fit to polynomial
  poly_param_fit<breit_wigner_KLOE> fit(rho);
  double starting_values[2] = {0., 0.,};
  fit.set_params(2, starting_values);
  fit.fit_params();

  // Print a Dalitz plot of the resulting fit
  test_poly fit_poly(1., 130., 20.);
  fit_poly.set_Mass(1.02);

  dalitz<test_poly> plot(fit_poly);
  plot.plot("POLY_dalitz_plot.dat");

  dalitz<breit_wigner_KLOE> plot_bw(rho);
  plot_bw.plot("BW_dalitz_plot.dat");

  return 0;
};
