#include "kt_amplitude.hpp"
#include "decay_kinematics.hpp"
#include "poly_exp.hpp"
#include "dalitz_fit.hpp"

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
  vector_meson.set_decayMass(1.02);
  vector_meson.set_decayIsospin(0);
  vector_meson.set_decayParticle("phi");

  // Options parameters for the KT euqations
  kt_options options;
  options.max_spin = 1;
  options.max_iters = 1;
  options.use_conformal = false;
  options.add_subtraction(1, 1, 1, 1);

  kt_amplitude kt_pwave(options, vector_meson);
  kt_pwave.iterate();
  kt_pwave.plot_inhomogeneity(0, 1);
  kt_pwave.sum_rule();
  //
  // kt_pwave.plot_isobar(0);
  // kt_pwave.plot_fundamentalSolutions(0);

  // // Start with the BESIII results
  // double pfit_1[2] = {1., 132.1};
  // poly_exp fit_1(2, pfit_1, vector_meson);
  // dalitz d_plot(&fit_1);
  //
  // // FIT 1
  // dalitz_fit fitter_1(&fit_1, &kt_pwave);
  // fitter_1.extract_params(2 * options.max_subs); // 2 real parameters = abs and arg of b
  // kt_pwave.print_params();
  // kt_pwave.normalize(7.56);
  // kt_pwave.plot_isobar(0);



  // // FIT 2
  // double pfit_2[3] = {1., 120.2, 29.5};
  // poly_exp fit_2(3, pfit_2, vector_meson);
  // // fit_2.normalize(7.56);
  //
  // dalitz_fit fitter_2(&fit_2, &kt_pwave);
  // fitter_2.extract_params(2 * options.max_subs);
  // // fitter_2.plot_deviation();
  // kt_pwave.print_params();
  // kt_pwave.plot_isobar(0);
  //
  // // FIT 3
  // double pfit_3[4] = {1., 111., 25., 22.};
  // poly_exp fit_3(4, pfit_3, vector_meson);
  // // fit_3.normalize(7.56);
  //
  // dalitz_fit fitter_3(&fit_3, &kt_pwave);
  // fitter_3.extract_params(2 * options.max_subs + 1);
  // kt_pwave.print_params();

  cout << endl;
  cout << "End." << endl;
  cout << "-----------------------------------------------------------" << endl;

  return 1.;
};
