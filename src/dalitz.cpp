// General routines for Dalitz plot generation
//
// Dependencies: constants.pp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dalitz.hpp"

//-----------------------------------------------------------------------------

// Doubly Differential Decay Width
template <class T>
double dalitz<T>::d2Gamma(double s, double t)
{
  complex<double> f;
  double temp1;
  temp1 = 96. * pow(2.* M_PI * amp.mDec, 3);
  f = dalitz<T>::amp(s,t);
  return abs(f * f) * amp.Kibble(s,t) / temp1;
};

//-----------------------------------------------------------------------------
// Quick function to print out a txt file with the Dalitz plot.
template <class T>
void dalitz<T>::plot(const char * n)
{
  std::ofstream output;
  output.open(n);

  double s_step = (amp.smax() - amp.smin() - offset)/100.;
  double t_step;

  double si, tij, gam;

  for (int i = 0; i < 100; i++)
  {
    si = amp.smin() + offset + double(i) * s_step;
    for(int j = 0; j < 100; j++)
    {
      t_step = (amp.tmax(si) - offset - amp.tmin(si))/100.;
      tij = amp.tmin(si) + offset  + double(j) * t_step;

      gam = d2Gamma(si, tij) ;

      output << std::left << setw(15) << amp.x(si, tij) << setw(15) << amp.y(si, tij) << setw << gam << endl;
    }
  }
output.close();
};
