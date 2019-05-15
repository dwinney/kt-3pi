#include "poly_param_fit.cpp"
#include "dalitz.cpp"
#include "breit_wigner.hpp"
#include "test_poly.hpp"

#include <iostream>
#include <complex>


using std::complex;
using std::cout;
using std::endl;

int main()
{
  double starting_values[3] = {0., 0., 0.};

  test_poly only_alpha(120., 29.5, 30.);
  only_alpha.set_Mass(1.02);

  dalitz<test_poly> test(only_alpha);
  test.plot("./BW_dalitz_plot.dat");
  //
  poly_param_fit<test_poly> testfit(only_alpha);
  testfit.set_params(2, starting_values);
  testfit.fit_params();
  // testfit.print_params();

  return 0;
};
