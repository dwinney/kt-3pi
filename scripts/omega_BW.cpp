#include "param_fit.cpp"
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
  breit_wigner rho(.770, .150, "rho_BW");
  rho.set_Mass(.78, "omega");

  dalitz<breit_wigner> plot_bw(rho);
  plot_bw.plot("BW_dalitz_plot.dat");

  // Fit to polynomial
  param_fit<breit_wigner> fit(rho);
  fit.extract_params(1);
  fit.extract_params(2);
  fit.extract_params(3);

  return 0;
};
