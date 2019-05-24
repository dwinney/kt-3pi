// Class definitions for KT formalism, including the rescattering integral and KT amplitude
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

using std::vector;
using std::complex;
using std::cout;


//-----------------------------------------------------------------------------

template <class T>
class kt_amplitude
//-----------------------------------------------------------------------------
{
protected:
T ini_lineshape;

//-----------------------------------------------------------------------------
public:
  kt_amplitude(T& my_lineshape) : ini_lineshape(my_lineshape){};

  
//-----------------------------------------------------------------------------
};
#endif
