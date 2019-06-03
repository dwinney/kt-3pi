// Class definitions for KT formalism, including the rescattering integral and KT amplitude
//
// Dependencies: amp.hpp
//
// Author:       Daniel Winney (2019)
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

// The kt class combines isobars together.
//-----------------------------------------------------------------------------
class kt : public isobar
//-----------------------------------------------------------------------------
{
protected:
vector<isobar> partial_waves;

//-----------------------------------------------------------------------------
public:
kt_amplitude() : ini_lineshape(my_lineshape){};


//-----------------------------------------------------------------------------
};
#endif
