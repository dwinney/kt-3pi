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

  // kt_pwave.print(0, 0);
  // kt_pwave.print(0, 1);
  kt_pwave.print_iteration(3, 0);
  kt_pwave.print_iteration(3, 1);

  // Start with the BESIII results
  double pfit1[2] = {1., 132.1};
  double pfit2[3] = {1., 120.2, 29.5};
  double pfit3[4] = {1., 111, 25, 22};

  poly_exp fit1(2, pfit1, vector_meson);
  poly_exp fit2(3, pfit2, vector_meson);
  poly_exp fit3(4, pfit3, vector_meson);

  // Now we fit!!
  dalitz_fit<poly_exp, isobar> fitter1(&fit1, &kt_pwave);
  fitter1.extract_params(options.max_subs + 1); // Remember there is a zeroth coefficient too!!!
  kt_pwave.print_params();
  kt_pwave.print();

  // dalitz_fit<poly_exp, isobar> fitter2(&fit2, &kt_pwave);
  // fitter2.extract_params(options.max_subs + 1);
  // kt_pwave.print_params();
  //
  // dalitz_fit<poly_exp, isobar> fitter3(&fit3, &kt_pwave);
  // fitter3.extract_params(options.max_subs + 1);
  // kt_pwave.print_params();

  return 1.;
};
