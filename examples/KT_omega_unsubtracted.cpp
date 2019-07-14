#include "kt_amplitude.hpp"
#include "dalitz_fit.cpp"
#include "dalitz.cpp"
#include "poly_exp.hpp"
#include "breit_wigner.hpp"

#include <iostream>
#include <complex>

using std::complex;
using std::cout;
using std::endl;

int main()
{
  cout << endl;

  // Set up the decay kinematics for the amplitude
  decay_kinematics vector_meson;
    vector_meson.set_decayJPC(1, -1, -1);
    vector_meson.set_decayMass(.780);
    vector_meson.set_decayIsospin(0);
    vector_meson.set_decayParticle("omega");

    cout << vector_meson.pseudo_threshold() << endl;

  // Options parameters for the KT euqations
  kt_options options;
  options.max_iters = 3;
  options.max_subs = 0;
  options.use_conformal = true;

  // Omega dominated by single isobar, j = 1, I = 1
  // in general this would be a KT object that is
  // defined with a Jmax and helicity / isospin dependence.
  isobar kt_pwave(1, 1, options, vector_meson);
  kt_pwave.iterate();
  kt_pwave.print_iteration(0, 0);
  kt_pwave.print_iteration(3, 0); // No rescattering

  cout << endl;
  cout << "Extracting Dalitz Plot Parameters..." << endl;
  cout << endl;

  poly_exp fit_results(vector_meson);

  dalitz_fit<isobar,poly_exp> fitter1(&kt_pwave, &fit_results);

  fitter1.extract_params(2);
  fit_results.print_params();

  fitter1.extract_params(3);
  fit_results.print_params();

  fitter1.extract_params(4);
  fit_results.print_params();

  // cout << "Extracting effective rho mass and width from a Breit-Wigner..." << endl;
  // cout << endl;
  //
  // breit_wigner rho_eff(vector_meson);
  //
  // dalitz_fit<isobar,breit_wigner> fitter2(&kt_pwave, &rho_eff);
  // fitter2.extract_params(2);
  // rho_eff.print_params();

  return 1.;
};
