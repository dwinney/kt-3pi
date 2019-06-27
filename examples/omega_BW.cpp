#include "dalitz_fit.cpp"
#include "dalitz.cpp"
#include "poly_exp.hpp"
#include "breit_wigner.hpp"
#include "omnes.hpp"

#include <iostream>
#include <complex>

using std::complex;
using std::cout;
using std::endl;

int main()
{
  // Set up the decay kinematics for the amplitude
  decay_kinematics vector_meson;
    vector_meson.set_decayJPC(1, 1, 1);
    vector_meson.set_decayMass(.780);
    vector_meson.set_decayIsospin(0);
    vector_meson.set_decayParticle("Omega");

  //----------------------------------------------------------------------------
  // Theoretical amplitude is a Breit-Wigner for the rho meson
  breit_wigner rho(.770, .150, vector_meson, "BW_physical_rho");
  rho.plot();

  dalitz<breit_wigner> plot_bw(&rho);
  plot_bw.plot();

  // Fit to polynomial expansion.
  poly_exp fit1(vector_meson); // Create empty poly_exp
  dalitz_fit<breit_wigner, poly_exp> fitter(&rho, &fit1); // Pass them both to the fitter

  fitter.extract_params(2); // Just alpha (normalization also a fit parameter)
  fit1.print_params();

  fitter.extract_params(3); // Both alpha and beta
  fit1.print_params();

  // Additionally we can start with the best fit BESIII values and calculate
  // the mass and width of and effective (rescattered) rho BW amplitude.
  double BES_params[3] = {4.78, 120.2, 29.5};
  poly_exp BES_fit(3, BES_params, vector_meson);

  breit_wigner rho_eff(vector_meson); //empty breit_wigner

  dalitz_fit<poly_exp, breit_wigner> fitter2(&BES_fit, &rho_eff);
  fitter2.extract_params(2);
  rho_eff.print_params();
  rho_eff.plot();

  return 0;
};
