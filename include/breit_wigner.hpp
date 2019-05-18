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

#include <complex>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

using std::setw;
using std::string;

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
  complex<double> f(double s);

// ---------------------------------------------------------------------------
};

// ---------------------------------------------------------------------------
// Reltivistic Breit-Wigner from the KLOE analysis [hep-ex/0204013]
// with finite-width corrections to VMD model of Rho decay
// ---------------------------------------------------------------------------
class breit_wigner: public amplitude
{
protected:
  string name;
  double res_mass, res_width;
  double s_res(){ return res_mass * res_mass;};
// ---------------------------------------------------------------------------
public:
  breit_wigner(double mass, double width, const char * n = "")
    : res_mass(mass), res_width(width), name(n)
  {
    cout << endl;
    if (name != "")
        {
      cout << name + ": ";
    }
    cout << "Breit-Wigner with M = " << mass << " and Gamma_0 = " << width << " created. \n";
  };

  complex<double> operator ()(double s, double t);
  complex<double> F(double x);
  complex<double> width(double s);
  complex<double> mom_pi(double s);

  void plot();
// ---------------------------------------------------------------------------
};
// ---------------------------------------------------------------------------

#endif