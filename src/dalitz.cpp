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
  complex<double> F = dalitz<T>::amp(s,t);
  double result = abs(amp.Kibble(s,t) * F * F);
  return result / normalization;
};

//-----------------------------------------------------------------------------
// Print out a txt file with the Dalitz plot and plot with ROOT
template <class T>
void dalitz<T>::plot(const char * filename)
{
  std::string name = filename;
  cout << "Plotting Dalitz plot... \n";

  std::ofstream output;
  name += ".dat";
  output.open(name.c_str());

  double s_step = (amp.smax() - amp.smin() - offset)/100.;
  double t_step;

  double si, tij, gam;

  for (int i = 0; i < 100; i++)
  {
    si = amp.smin() + offset + double(i) * s_step;
    for(int j = 0; j < 100; j++)
    {
      t_step = (amp.tmax(si) - offset - amp.tmin(si)) / 100.;
      tij = amp.tmin(si) + offset  + double(j) * t_step;

      gam = d2Gamma(si, tij);

      output << std::left << setw(15) << amp.x(si, tij) << setw(15) << amp.y(si, tij) << setw(15) << gam << endl;
    }
  }
output.close();

cout << "Output to : " << name << endl;
cout << endl;

TCanvas *c = new TCanvas("c", "c");
TGraph2D *g = new TGraph2D(name.c_str());
g->Draw("colz");
gStyle->SetPalette(55);

c->Modified();
name.erase(name.end() - 4, name.end());
name += ".pdf";
c->Print(name.c_str());
c->Destructor();
};
