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
  rho.plot();

  dalitz<breit_wigner> plot_bw(rho);
  plot_bw.plot();

  // Fit to polynomial expansion.
  dalitz_fit<breit_wigner, poly_exp> fitter(rho);
  poly_exp fit1 = fitter.extract_params(2); // Just alpha (normalization also a fit parameter)
  fit1.print_params();

  poly_exp fit2 = fitter.extract_params(3); // Both alpha and beta
  fit2.print_params();

  cout << "--------------------------------------------------------------------- \n";
  // Additionally we can start with the best fit BESIII values and calculate
  // the mass and width of and effective (rescattered) rho BW amplitude.
  double BES_params[3] = {4.78, 120.2, 29.5};
  poly_exp BES_fit(3, BES_params);
  BES_fit.set_decayMass(.780);

  dalitz_fit<poly_exp, breit_wigner> fitter2(BES_fit);
  breit_wigner rho_eff = fitter2.extract_params(2);
  rho_eff.print_params();
  rho_eff.plot();

    cout << "--------------------------------------------------------------------- \n";
  //----------------------------------------------------------------------------
  // We can do this with an arbitrary input function.
  // Here we compare with the pure Omnes function, extracting Dalitz plot parameters.
  omnes pwave(1., 1., "Omnes_pwave"); // isospin 1, spin 1
  pwave.set_decayMass(.780);

  dalitz<omnes> plot_omnes(pwave);
  plot_omnes.plot();

  dalitz_fit<omnes, poly_exp> fitter3(pwave);
  poly_exp Omnes_fit = fitter3.extract_params(2);
  Omnes_fit.print_params();
  //----------------------------------------------------------------------------
  return 0;
};
