#include "kt_amplitude.hpp"

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

  // Omega dominated by single isobar, in general this would be a KT object that is
  // defined with a Jmax and helicity / isospin dependence.
  isobar pwave(1, 1, vector_meson);
  pwave.iterate(3); // Calculate three iterations of the KT equation

  // Print the results compared to the no-3body-effects case.
  pwave.print(0); // No rescattering
  pwave.print(1);
  pwave.print(3); // 3rd iterations

return 1.;
};