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

  // Similar to the unsubtracted case we start by calculating the dispersion relations
  // Options parameters for the KT euqations
  kt_options options;
  options.max_iters = 3;
  options.max_subs = 1;
  options.use_conformal = true;

  isobar kt_pwave(1, 1, options, vector_meson);
  kt_pwave.iterate(); // Calculate three iterations of the KT equation

  // kt_pwave.print_iteration(3, 0);
  // kt_pwave.print_iteration(3, 1);

  // Start with the BESIII results
  double pfit_1[2] = {1., 132.1};
  double perr_1[2] = {0., 6.7};
  poly_exp fit_1(2, pfit_1, vector_meson);
  fit_1.set_errors(2, perr_1);

  // FIT 1
  dalitz_fit<poly_exp, isobar> fitter_1(&fit_1, &kt_pwave);
  fitter_1.extract_params(options.max_subs + 1); // Remember there is a zeroth coefficient too!!!
  kt_pwave.print_params();
  kt_pwave.print();

  // FIT 2
  double pfit_2[3] = {1., 120.2, 29.5};
  double perr_2[3] = {0., 7.1, 8.0};
  poly_exp fit_2(3, pfit_2, vector_meson);
  fit_2.set_errors(3, perr_2);

  dalitz_fit<poly_exp, isobar> fitter_2(&fit_2, &kt_pwave);
  fitter_2.extract_params(options.max_subs + 1);
  kt_pwave.print_params();

  // FIT 3
  double pfit_3[4] = {1., 111., 25., 22.};
  double perr_3[4] = {0., 18., 10., 29.};
  poly_exp fit_3(4, pfit_3, vector_meson);
  fit_3.set_errors(4, perr_3);

  dalitz_fit<poly_exp, isobar> fitter_3(&fit_3, &kt_pwave);
  fitter_3.extract_params(options.max_subs + 1);
  kt_pwave.print_params();

  return 1.;
};
