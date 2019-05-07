// Classes for General Amplitude Stuff including kinematic functions such as momenta.
//
// Dependencies: amp.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amp.hpp"

// Mandelstam u in terms of the other two variables
double amplitude::u(double s, double t)
{
 return 3.*mPi*mPi + Mass*Mass - s - t;
}

// Lorentz Covariant Kibble function, only defined in the physical reqgions, i.e. it can only be real
double amplitude::Kibble(double s, double t)
{
    double u = 3.*mPi*mPi + Mass*Mass - s - t;
    return s * t * u - mPi*mPi* std::pow(Mass*Mass - mPi*mPi, 2);
};
