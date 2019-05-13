#include "poly_param_fit.hpp"
#include "breit_wigner.hpp"
#include <iostream>
#include <complex>
#include <functional>

using std::complex;
using std::cout;
using std::endl;

int main()
{
  // using namespace std::placeholders;
  breit_wigner_simple rho(.770, .150);

  // std::function<complex<double>(double, double)> rho_BW = std::bind(&breit_wigner_simple::operator(), rho, _1, _2);

  // test.set_Mass(.782);
  // test.set_params(1., a);
  // test.fit_params();
  // test.print_params();

  return 0;
};
