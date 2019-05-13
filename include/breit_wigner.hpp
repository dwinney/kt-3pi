// Class for different Breit-Wigner parameterizations of amplitudes
//
// Dependencies: amp.cpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _BW_
#define _BW_

#include "amp.hpp"

// ---------------------------------------------------------------------------
// Relativisitic Breit-Wigner with a constant imaginary part.
// ---------------------------------------------------------------------------
class breit_wigner_simple : public amplitude
{
protected:
  double res_mass, res_width;
  double s_res(){ return res_mass * res_mass;};
// ---------------------------------------------------------------------------
public:
  breit_wigner_simple(double mass, double width)
  {
  res_mass = mass;
  res_width = width;
  };

  complex<double> operator ()(double s, double t);
// ---------------------------------------------------------------------------
};

// ---------------------------------------------------------------------------
// Reltivistic Breit-Wigner from the KLOE analysis [hep-ex/0204013]
// with an Energy-dependent width.
// ---------------------------------------------------------------------------
class breit_wigner_KLOE : public amplitude
{
protected:
  double res_mass, res_width;
  double s_res(){ return res_mass * res_mass;};
// ---------------------------------------------------------------------------
public:
  breit_wigner_KLOE(double mass, double width)
  {
  res_mass = mass;
  res_width = width;
  };

  complex<double> operator ()(double s, double t);
  double width(double s);
// ---------------------------------------------------------------------------
};
// ---------------------------------------------------------------------------

#endif
