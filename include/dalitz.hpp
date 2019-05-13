// General purpose classes for Dalitz plot generation.
// for fitting to polynomial expansion see poly_param_fit.hpp
//
// Dependencies: amp.hpp, aux_math.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _DALITZ_
#define _DALITZ_

#include "amp.hpp"
#include <cmath>

//-----------------------------------------------------------------------------
//  Define a Dalitz plot object with a pointer to the amplitude.
//  The amp should be a complex<double> and a function of 2 doubles, Mandelstam s and t.
//
//  complex<double> (*ptr_amp) (double, double);
//  ptr_func = my_amp;
//
//  dalitz my_dalitz(ptr_amp);
//-----------------------------------------------------------------------------

template <class T>
class dalitz : public amplitude
{
protected:
T amp;

// Center of the Dalitz plot in s and t
double s_c()
  {
    return  (mDec*mDec + 3. * mPi*mPi) / 3.;
  };
double t_c()
{
    return (3.*mPi*mPi + mDec*mDec - s_c())/2.;
};

//Temporary variables
double temp1, temp2, temp3;

//-----------------------------------------------------------------------------
public:
  dalitz()
  {
    cout << " Need pointer to amplitude function in dalitz object declaration. Quitting... \n";
    std::exit(1);
  };
  dalitz(T& my_amp) : amp(my_amp){};
  complex<double> operator()(double s, double t)
  {
    return amp(s,t);
  };
//-----------------------------------------------------------------------------
  // Lorentz Invariant dimensionless parameters
  double x(double s, double t);
  double y(double s, double t);

  double x_polar(double z, double theta);
  double y_polar(double z, double theta);

  double z(double s, double t);
  double theta(double s, double t);

  // Double differential cross section
  double d2Gamma(double s, double t);
};

//-----------------------------------------------------------------------------
#endif
