#include "dalitz_fit.cpp"
#include "dalitz.cpp"
#include "breit_wigner.hpp"
#include "omnes.hpp"
#include "kt_isobar.hpp"

#include <iostream>
#include <complex>

#include "Math/Interpolator.h"

using std::complex;
using std::cout;
using std::endl;

int main()
{
  vector<double> s, fx;
  breit_wigner rho(.770,.150);

  for (int i = 0; i < 200; i++)
  {
    double s_i = sthPi + double(i)* (elastic_cutoff - sthPi) / 200.;

    s.push_back(s_i);

    complex<double> f_rho = rho.F(s_i);
    fx.push_back( real(f_rho) );
  }

  ROOT::Math::Interpolator inter(s, fx, ROOT::Math::Interpolation::kCSPLINE) ;

  cout << std::left << std::setw(15) << real(rho.F(.666)) << std::setw(15) << inter.Eval(.666) << std::endl;

 return 1.;
};
