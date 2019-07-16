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
#include "dalitz.cpp"

#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TError.h>

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
  breit_wigner_simple(decay_kinematics dec)
  : kinematics(dec)
  {};

  breit_wigner_simple(double mass, double width, decay_kinematics dec)
  : kinematics(dec)
  {
  res_mass = mass;
  res_width = width;
  };

  decay_kinematics kinematics;

  complex<double> eval(double s, double t);
  complex<double> f(double s);
  double error_func(double s, double t)
  {
    return 1.;
  };

// ---------------------------------------------------------------------------
};

// ---------------------------------------------------------------------------
// Reltivistic Breit-Wigner from the KLOE analysis [hep-ex/0204013]
// with finite-width corrections to VMD model of Rho decay
// ---------------------------------------------------------------------------
class breit_wigner
{
protected:
  double normalization = 1.;
  double res_mass, res_width;
  double s_res(){ return res_mass * res_mass;};

  complex<double> width(double s);
  complex<double> mom_pi(double s);
// ---------------------------------------------------------------------------
public:
  breit_wigner(decay_kinematics dec)
  : kinematics(dec)
  {};

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
  complex<double> eval(double s, double t);
  double error_func(double s, double t)
  {
    return 1.;
  };

  void normalize(double gamma_exp);
  void set_params(int n, const double * par);
  void print_params();
  void print();
// ---------------------------------------------------------------------------
};
// ---------------------------------------------------------------------------

#endif
