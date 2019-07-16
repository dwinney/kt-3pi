#include "kt_amplitude.hpp"
#include "dalitz_fit.cpp"
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
  kt_pwave.print_iteration(options.max_iters, 0);

  cout << endl;
  cout << "Extracting Dalitz Plot Parameters..." << endl;
  cout << endl;

  poly_exp fit_results(vector_meson);

  dalitz_fit<isobar,poly_exp> fitter1(&kt_pwave, &fit_results);
  kt_pwave.normalize(7.56);
  kt_pwave.print();

  fitter1.extract_params(2);
  fit_results.print_params();
  //
  fitter1.extract_params(3);
  fit_results.print_params();

  fitter1.extract_params(4);
  fit_results.print_params();

  return 1.;
};
