// Class for different Breit-Wigner parameterizations of amplitudes
//
// Dependencies: decay_kinematics.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _BW_
#define _BW_

#include "decay_kinematics.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>

using std::setw;

// ---------------------------------------------------------------------------
// Relativisitic Breit-Wigner with a constant imaginary part.
class breit_wigner_simple
{
protected:
  double res_mass, res_width;
  double s_res(){ return res_mass * res_mass;};
// ---------------------------------------------------------------------------
public:
  breit_wigner_simple(double mass, double width, decay_kinematics dec)
  : kinematics(dec)
  {
  res_mass = mass;
  res_width = width;
  };

  decay_kinematics kinematics;

  complex<double> operator ()(double s, double t);
  complex<double> f(double s);

// ---------------------------------------------------------------------------
};

// ---------------------------------------------------------------------------
// Reltivistic Breit-Wigner from the KLOE analysis [hep-ex/0204013]
// with finite-width corrections to VMD model of Rho decay
// ---------------------------------------------------------------------------
class breit_wigner
{
protected:
  double res_mass, res_width;
  double s_res(){ return res_mass * res_mass;};

  complex<double> width(double s);
  complex<double> mom_pi(double s);
// ---------------------------------------------------------------------------
public:
  breit_wigner(){};
  breit_wigner(double mass, double width, decay_kinematics dec, const char * n = "")
    : res_mass(mass), res_width(width), kinematics(dec)
  {
    kinematics.set_ampName(n);
    cout << endl;
    if (kinematics.get_ampName() != "")
        {
      cout << kinematics.get_ampName()  + ": ";
    }
    cout << "Breit-Wigner with M = " << mass << " and Gamma_0 = " << width << " created. \n";
  };

  decay_kinematics kinematics;

  complex<double> F(double x);
  complex<double> operator ()(double s, double t);

  void set_params(int n, const double * par);
  void print_params();
  void plot();
// ---------------------------------------------------------------------------
};
// ---------------------------------------------------------------------------

#endif
