#include "poly_param_fit.cpp"
#include "poly_param_fit.hpp"
#include "dalitz.hpp"
#include "dalitz.cpp"
#include "breit_wigner.hpp"
#include <iostream>
#include <complex>

using std::complex;
using std::cout;
using std::endl;

int main()
{
  breit_wigner_KLOE rho(.770, .150);
  poly_param_fit<breit_wigner_KLOE> test(rho);
  test.set_Mass(.782);
  test.set_integration_points(100);
  double start[1] = {0.};
  test.set_params(1, start);
  test.fit_params();
  test.print_params();

  return 0;
};
