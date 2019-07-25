#include "decay_kinematics.hpp"
#include "poly_exp.hpp"
#include "breit_wigner.hpp"

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
    vector_meson.set_decayParticle("omega");

  //----------------------------------------------------------------------------
  // Theoretical amplitude is a Breit-Wigner for the rho meson
  breit_wigner rho(.770, .150, vector_meson, "BW_physical_rho");
  rho.normalize(7.56);
  rho.print();

  // Fit to polynomial expansion.
  poly_exp fit(vector_meson); // Create empty poly_exp
  dalitz_fit<breit_wigner, poly_exp> fitter(&rho, &fit); // Pass them both to the fitter

  fitter.extract_params(2); // Just alpha (normalization also a fit parameter)
  fit.print_params();
  fitter.plot();

  fitter.extract_params(3); // Both alpha and beta
  fit.print_params();

  fitter.extract_params(4); // Both alpha and beta
  fit.print_params();

  // // Additionally we can start with the best fit BESIII values and calculate
  // // the mass and width of and effective (rescattered) rho BW amplitude.
  // double BES_params[3] = {1., 120.2, 29.5};
  // poly_exp BES_fit(3, BES_params, vector_meson);
  //
  // dalitz<poly_exp> test(&BES_fit);
  // test.plot();
  //
  // breit_wigner rho_eff(vector_meson); //empty breit_wigner
  // dalitz_fit<poly_exp, breit_wigner> fitter2(&BES_fit, &rho_eff);
  // fitter2.extract_params(2);
  // rho_eff.normalize(7.56);
  // rho_eff.print_params();
  // rho_eff.print();

  return 0;
};
