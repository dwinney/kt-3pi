#include "kt_amplitude.hpp"

#include <iostream>
#include <complex>

#include "Math/Interpolator.h"

using std::complex;
using std::cout;
using std::endl;

int main()
{

decay_kinematics omega;
  omega.set_decayJPC(1, 1, 1);
  omega.set_decayMass(1.1);
  omega.set_decayIsospin(0);
  omega.set_decayParticle("Omega");
  
return 1.;
};
