#include "param_fit.cpp"
#include "dalitz.cpp"
#include "breit_wigner.hpp"
#include "omnes.hpp"

#include <iostream>
#include <complex>

using std::complex;
using std::cout;
using std::endl;

int main()
{
  //----------------------------------------------------------------------------
  // Theoretical amplitude is a Breit-Wigner for the rho meson
  breit_wigner rho(.770, .150, "BW_physical_rho");
  rho.set_decayMass(.780);

  dalitz<breit_wigner> plot_bw(rho);
  plot_bw.plot();

  // Fit to polynomial
  param_fit<breit_wigner> fit1(rho);
  fit1.extract_params(1); // Just alpha
  fit1.extract_params(2); // Both alpha and beta

  //----------------------------------------------------------------------------
  //Additionally we can compare with the pure Omnes function
  omnes pwave(1., 1., "Omnes_pwave"); // isospin 1, spin 1
  pwave.set_decayMass(.780);

  dalitz<omnes> plot_omnes(pwave);
  plot_omnes.plot();

  param_fit<omnes> fit2(pwave);
  fit2.extract_params(1);
  //----------------------------------------------------------------------------
  return 0;
};
