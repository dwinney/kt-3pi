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

#include "kt_isobar.hpp"

// The kt class combines isobars together.
//-----------------------------------------------------------------------------
class kt_amplitude
//-----------------------------------------------------------------------------
{
protected:
int max_iter;

vector<isobar> isobars;

//-----------------------------------------------------------------------------
public:
kt_amplitude(){};


//-----------------------------------------------------------------------------
};
#endif
