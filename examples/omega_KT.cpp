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
    vector_meson.set_decayParticle("Omega");

  // Omega dominated by single isobar, j = 1, I = 1
  // in general this would be a KT object that is
  // defined with a Jmax and helicity / isospin dependence.
  int n_subtractions = 0;
  isobar kt_pwave(1, 1, n_subtractions, vector_meson);
  kt_pwave.iterate(3); // Calculate three iterations of the KT equation

  kt_pwave.print(0, 0); // No rescattering
  kt_pwave.print(3, 0);

  // Extract dalitz plot parameters
  dalitz_fit<isobar,poly_exp> fitter(kt_pwave);
  poly_exp fit_results = fitter.extract_params(3);
  fit_results.print_params();

  // Compare with a Breit Wigner
  dalitz_fit<isobar,breit_wigner> fitter2(kt_pwave);
  breit_wigner BW_fit = fitter2.extract_params(2);
  BW_fit.print_params();

return 1.;
};
