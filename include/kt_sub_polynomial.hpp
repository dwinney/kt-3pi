
#ifndef _DISP_POLY_
#define _DISP_POLY_

#include "constants.hpp"

// ---------------------------------------------------------------------------
// Super simple class to output an nth order polynomial in a conformal variable
// ---------------------------------------------------------------------------
class subtraction_polynomial
{
protected:
  const complex<double> s_inelastic = xr;
  const complex<double> s_expand = 0.;
public:
  subtraction_polynomial(){};

  // ----------------------------------------------------------------------------
  // Conformal variable which maps the cut plane in complex s to the unit disk
  complex<double> conformal(complex<double> s, int ieps)
  {
    complex<double> numerator, denominator;

    numerator = sqrt(s_inelastic - s_expand) - sqrt(s_inelastic - s + xi * double(ieps) * EPS);
    denominator = sqrt(s_inelastic - s_expand) + sqrt(s_inelastic - s + xi * double(ieps) * EPS);

    return numerator / denominator;
  };

  // ----------------------------------------------------------------------------
  // Outputs a polynomial of order n with unit coefficients in the above conformal variable
  complex<double> operator() (int n, complex<double> s, int ieps)
  {
    return pow(conformal(s, ieps), double(n));
  };

};

#endif
