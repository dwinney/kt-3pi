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
#include "dalitz.hpp"

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
// Reltivistic Breit-Wigner from the KLOE analysis [hep-ex/0204013]
// with finite-width corrections to VMD model of Rho decay
// ---------------------------------------------------------------------------

class breit_wigner : public amplitude
{
protected:
  double normalization = 1.;
  double res_mass, res_width;
  double s_res(){ return res_mass * res_mass;};

  complex<double> width(double s);
  complex<double> mom_pi(double s);
  complex<double> F(double x);
  
// ---------------------------------------------------------------------------
public:
  breit_wigner(decay_kinematics dec)
  : amplitude(dec)
  {};

  breit_wigner(double mass, double width, decay_kinematics dec, const char * n = "")
    : res_mass(mass), res_width(width), amplitude(dec)
  {
    dec.set_ampName(n);
    cout << endl;
    if (dec.get_ampName() != "")
        {
      cout << dec.get_ampName()  + ": ";
    }
    cout << "Breit-Wigner with M = " << mass << " and Gamma_0 = " << width << " created. \n";
  };

  complex<double> eval(double s, double t);

  void normalize(double gamma_exp);
  void set_params(int n, const double * par);
  void print_params();
  void print();

};
// ---------------------------------------------------------------------------

#endif
