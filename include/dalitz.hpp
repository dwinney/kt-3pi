// Classes for Dalitz plot generation and fitting to polynomial expansion with dalitz plot parameters.
//
//
// Dependencies: amp.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _DALITZ_
#define _DALITZ_

#include "amp.hpp"
#include <iostream>
#include <complex>
#include <vector>
using std::vector;
using std::complex;
using std::cout;

//-----------------------------------------------------------------------------
//  Define a Dalitz plot object with a pointer to the amplitude.
//  The amp should be a complex<double> and a function of 2 doubles, Mandelstam s and t.
//
//  complex<double> (*ptr_amp) (double, double);
//  ptr_func = my_amp;
//
//  dalitz my_dalitz(ptr_amp);
//-----------------------------------------------------------------------------

class dalitz : amplitude
{
protected:
complex<double> (*amp) (double, double);

public:
  dalitz()
  {
    cout <<
    " Need pointer to amplitude function in datlitz object declaration. Quitting..."
     << exit(1);
  }
  dalitz(complex<double> (*my_amp) (double, double))
  {
    amp = my_amp;
  }

  double kibble(double s, double t);

};

//-----------------------------------------------------------------------------
#endif
