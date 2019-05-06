// Class definitions for KT formalism, including the rescattering integral and KT amplitude
//
//
// Dependencies: amp.hpp
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _KT_
#define _KT_

#include <iostream>
#include <complex>
#include <vector>
#include "amp.hpp"
using std::vector;
using std::complex;
using std::cout;


//-----------------------------------------------------------------------------
//  Define a KT amplitude with a pointer to  initial funtion which will into the integral (i.e. a lineshape with no rescattering)
//  Function needs two arguments, the energy as a double and epsilon (+1 for above cut, -1 for below cut, 0 on the cut) as an int.
//
//  complex<double> (*ptr_func) (double, int);
//  ptr_func = my_function;
//
//  kt_amplitude my_kt(ptr_func);
//-----------------------------------------------------------------------------
class kt_amplitude : amplitude
//-----------------------------------------------------------------------------
{
protected:
complex<double> (*ini_lineshape) (double, int);
int N_points, N_iter; //number of iterations in rescattering equation
vector<complex<double>> interp_val;
//-----------------------------------------------------------------------------
public:
  kt_amplitude(complex<double> (*my_lineshape)(double, int), int n = 200)
  {
    ini_lineshape = my_lineshape;
    N_points = n;
  };

  void start();

  complex<double> eval(double s, int epsilon)
  {
    return ini_lineshape(s, epsilon);
  };
//-----------------------------------------------------------------------------
};
#endif
