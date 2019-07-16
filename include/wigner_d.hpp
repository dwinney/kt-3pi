// Functions related to the Wigner little-d functions and the half-angle factors
// which remove the kinematic singularities in the cross-variable
//
// Dependencies: None
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#ifndef _D_FUNC_
#define _D_FUNC_

#include <TMath.h>
#include <Math/SpecFuncMathMore.h>

// -----------------------------------------------------------------------------
// Outputs the Wigner little-d function appropriate for the decay into three scalar,
// d^j_{lambda, 0}(z)
double d_func(int j, int l, double z);

// -----------------------------------------------------------------------------
// Outputs d_hat the kinematic-singularity-free d_function
double d_hat(int j, int l, double z);

#endif
