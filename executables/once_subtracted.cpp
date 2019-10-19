#include "kt_amplitude.hpp"
#include "decay_kinematics.hpp"
#include "dalitz.hpp"
#include "dalitz_fit.hpp"
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
  cout << "-----------------------------------------------------------" << endl;

  // Set up the decay kinematics for the amplitude
  decay_kinematics vector_meson;
    vector_meson.set_decayJPC(1, -1, -1);
    vector_meson.set_decayMass(0.78);
    vector_meson.set_decayParticle("omega");

  // Options parameters for the KT equations
  kt_options options;
  options.max_iters = 1;
  options.max_subs = 0;
  options.max_spin = 1;
  options.interp_cutoff = 1.3;
  options.use_conformal = false;
  // options.test_angular = true;

  kt_amplitude kt_pwave(options, vector_meson);
  kt_pwave.iterate();

  kt_pwave.normalize(8.49);
  kt_pwave.plot_isobar(0);

  dalitz d_plot(&kt_pwave);
  d_plot.plot("KSF normalized");
  d_plot.plot("normalized");

  cout << "-----------------------------------------------------------" << endl;

  cout << "Extracting Dalitz Plot Parameters..." << endl;
  cout << endl;

  poly_exp fit_results(vector_meson);

  dalitz_fit fitter(&kt_pwave, &fit_results);

  fitter.extract_params(2);
  fit_results.print_params();
  fitter.plot_deviation();

  fitter.extract_params(3);
  fit_results.print_params();
  fitter.plot_deviation();

  fitter.extract_params(4);
  fit_results.print_params();
  fitter.plot_deviation();

  cout << "End." << endl;
  cout << "-----------------------------------------------------------" << endl;

  return 1.;
};
