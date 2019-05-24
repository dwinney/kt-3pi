#include "param_fit.cpp"
#include "dalitz.cpp"
#include "breit_wigner.hpp"
#include "test_poly.hpp"
#include "omnes.hpp"

#include <iostream>
#include <complex>

using std::complex;
using std::cout;
using std::endl;

int main()
{
 // Theoretical amplitude is a Breit-Wigner for the rho
 // in each channel
  breit_wigner rho(.770, .150);
  rho.set_ampName("BW_physical_rho");
  rho.set_decayMass(.78);


  dalitz<breit_wigner> plot_bw(rho);
  plot_bw.plot();

  // Fit to polynomial
  param_fit<breit_wigner> fit(rho);
  fit.extract_params(1);
  fit.extract_params(2);

  //Additionally we can compare with the Omnes solution
  omnes pwave(1., 1.);
  pwave.set_decayMass(.78);
  pwave.set_ampName("Omnes_pwave");

  dalitz<omnes> plot_omnes(pwave);
  plot_omnes.plot();

  param_fit<omnes> fit2(pwave);
  fit2.extract_params(1);

  return 0;
};
