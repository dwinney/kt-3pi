#include "kt_amplitude.hpp"
#include "decay_kinematics.hpp"
#include "poly_exp.hpp"
#include "breit_wigner.hpp"

#include <iostream>
#include <complex>
#include <array>

using std::complex;
using std::cout;
using std::endl;

int main()
{
  cout << endl;
  cout << "-----------------------------------------------------------" << endl;

  // Set up the decay kinematics for the amplitude
  decay_kinematics vector_meson;
  vector_meson.set_decayJPC(1, -1, -1);
  vector_meson.set_decayMass(.780);
  vector_meson.set_decayIsospin(0);
  vector_meson.set_decayParticle("omega");

  // Options parameters for the KT euqations
  kt_options options;
  options.max_spin = 1;
  options.max_iters = 1;
  options.use_conformal = false;
  options.add_subtraction(1, 1, 1, 1);

  kt_amplitude kt_pwave(options, vector_meson);
  kt_pwave.iterate();

  // kt_pwave.plot_iteration(1, 0, 0);
  // kt_pwave.plot_iteration(1, 0, 1);

  // // Start with the BESIII results
  // double pfit_1[2] = {1., 132.1};
  // double perr_1[2] = {0., 6.7};
  // poly_exp fit_1(2, pfit_1, vector_meson);
  // fit_1.normalize(7.56);
  // fit_1.set_errors(2, perr_1);
  //
  // // FIT 1
  // dalitz_fit<poly_exp, isobar> fitter_1(&fit_1, &kt_pwave);
  // fitter_1.extract_params(2 * options.max_subs + 1); // 2 real parameters = abs and arg of b
  // kt_pwave.print_params();
  //
  // // FIT 2
  // double pfit_2[3] = {1., 120.2, 29.5};
  // double perr_2[3] = {0., 7.1, 8.0};
  // poly_exp fit_2(3, pfit_2, vector_meson);
  // fit_2.normalize(7.56);
  // fit_2.set_errors(3, perr_2);
  //
  // dalitz_fit<poly_exp, isobar> fitter_2(&fit_2, &kt_pwave);
  // fitter_2.extract_params(2 * options.max_subs + 1);
  // kt_pwave.print_params();
  // kt_pwave.print();
  //
  // // FIT 3
  // double pfit_3[4] = {1., 111., 25., 22.};
  // double perr_3[4] = {0., 18., 10., 29.};
  // poly_exp fit_3(4, pfit_3, vector_meson);
  // fit_3.normalize(7.56);
  // fit_3.set_errors(4, perr_3);
  //
  // dalitz_fit<poly_exp, isobar> fitter_3(&fit_3, &kt_pwave);
  // fitter_3.extract_params(2 * options.max_subs + 1);
  // kt_pwave.print_params();

  cout << endl;
  cout << "End." << endl;
  cout << "-----------------------------------------------------------" << endl;

  return 1.;
};
