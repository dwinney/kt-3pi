#include "dalitz_fit.cpp"
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
  dalitz_fit<breit_wigner, poly_exp> fitter(rho);
  poly_exp fit1 = fitter.extract_params(2); // Just alpha (normalization also a fit parameter)
  fit1.print_params();

  poly_exp fit2 = fitter.extract_params(3); // Both alpha and beta
  fit2.print_params();


  //----------------------------------------------------------------------------
  // //Additionally we can compare with the pure Omnes function
  // omnes pwave(1., 1., "Omnes_pwave"); // isospin 1, spin 1
  // pwave.set_decayMass(.780);
  //
  // dalitz<omnes> plot_omnes(pwave);
  // plot_omnes.plot();
  //
  // param_fit<omnes> fit2(pwave);
  // fit2.extract_params(1);
  //----------------------------------------------------------------------------
  return 0;
};
