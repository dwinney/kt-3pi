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

  // Set up the decay kinematics for the amplitude
  decay_kinematics vector;
    vector.set_decayJPC(1, 1, 1);
    vector.set_decayMass(.780);
    vector.set_decayIsospin(0);
    vector.set_decayParticle("Omega");

  //----------------------------------------------------------------------------
  // Theoretical amplitude is a Breit-Wigner for the rho meson
  breit_wigner rho(.770, .150, vector, "BW_physical_rho");
  rho.plot();

  dalitz<breit_wigner> plot_bw(rho);
  plot_bw.plot();

  // Fit to polynomial expansion.
  dalitz_fit<breit_wigner, poly_exp> fitter(rho);
  poly_exp fit1 = fitter.extract_params(2); // Just alpha (normalization also a fit parameter)
  fit1.print_params();

  poly_exp fit2 = fitter.extract_params(3); // Both alpha and beta
  fit2.print_params();

  cout << "-------------------------------------------------------------------- \n";
  cout << endl;

  // Additionally we can start with the best fit BESIII values and calculate
  // the mass and width of and effective (rescattered) rho BW amplitude.
  double BES_params[3] = {4.78, 120.2, 29.5};
  poly_exp BES_fit(3, BES_params, vector);

  dalitz_fit<poly_exp, breit_wigner> fitter2(BES_fit);
  breit_wigner rho_eff = fitter2.extract_params(2);
  rho_eff.print_params();
  rho_eff.plot();

  cout << "-------------------------------------------------------------------- \n";
  cout << endl;

  return 0;
};
