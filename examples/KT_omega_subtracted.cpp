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

    // Similar to the unsubtracted case we start by calculating the dispersion relations
    int n_subtractions = 1;
    isobar kt_pwave(1, 1, n_subtractions, vector_meson);
    kt_pwave.iterate(3); // Calculate three iterations of the KT equation

    kt_pwave.print(3, 0);
    kt_pwave.print(3, 1);

    // Start with the BESIII results
    double BES_params[3] = {1., 120.2, 29.5};
    poly_exp BES_fit(3, BES_params, vector_meson);

    // Now we fit!!
    dalitz_fit<poly_exp, isobar> fitter(&BES_fit, &kt_pwave);
    fitter.extract_params(n_subtractions + 1); // Remember there is a zeroth coefficient too!!!

    kt_pwave.print_params();

  return 1.;
};
