#include "kt_amplitude.hpp"

#include <iostream>
#include <complex>

using std::complex;
using std::cout;
using std::endl;

int main()
{

decay_kinematics vector_meson;
  vector_meson.set_decayJPC(1, 1, 1);
  vector_meson.set_decayMass(0.780);
  vector_meson.set_decayIsospin(0);
  vector_meson.set_decayParticle("Omega");

isobar pwave(1, 1, vector_meson);
pwave.iterate(3);


return 1.;
};
