#include "poly_param_fit.hpp"
#include <iostream>
#include <complex>

using std::complex;
using std::cout;

complex<double> test_func(double s, double t)
{
  return xr;
};

int main()
{
  complex<double> (*ptr_amp) (double, double);
  ptr_amp = test_func;

  double a[2] = {0., 0.};

  poly_param_fit test(test_func);
  test.set_Mass(.782);
  // test.print_params();
  test.set_params(2, a);
  test.print_params();
  test.set_integration_points(10);
  cout << test.chi_squared() << "\n";
  return 0;
};
