// These objects include the isobars (i.e. the partial waves in each channel)
// they inherite from the amplitude class because they describe the scattering dynamics
// of a single spin and isospin-projected subchannel.
//
// Dependencies: amp.hpp, omnes.hpp, iteration.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "kt_isobar.hpp"

// ----------------------------------------------------------------------------
// Populate the zeroth iteration with values directly from call to Omnes function
void isobar::start()
{
  cout << endl;
  cout << "Storing initial Omnes amplitude... ";
  vector<double> s;

  // Need to store values of the Omnes amplitude around the unitarity cut for the
  // analytic continuation of the angular integral.
  vector< complex<double> > above_cut, below_cut;
  complex<double> ab, be;

  for (int i = 0; i < interpolation::N_interp; i++)
  {
    double s_i = sthPi + EPS + double(i) * (elastic_cutoff - sthPi - EPS) / double(interpolation::N_interp);
    s.push_back(s_i);

    omega.set_ieps(1); // above the unitarity cut (+ i epsilon)
    ab = omega.eval(s_i);
    above_cut.push_back(ab);

    omega.set_ieps(-1); // below the cut (- i epsilon)
    be = omega.eval(s_i);
    below_cut.push_back(be);
  }

  // If there are previously stored iterations for some reason, clear it when start() is called.
  if (iters.size() != 0)
  {
    iters.clear();
    cout << "isobar: Clearing previously stored Iterations for j = " << spin_proj;
    cout << " and I = " << iso_proj << ". \n";
  }

  iteration zeroth(0, omega, s, above_cut, below_cut);
  iters.push_back(zeroth);
  cout << "Done." << endl;
};

// ----------------------------------------------------------------------------
// Iterate through the KT equations n times, storing each iteration for comparison
void isobar::iterate(int n)
{
  start();

  for (int i = 1; i < n + 1; i++)
  {
  cout << "Calculating iteration (" << i << "/" << n << ")... ";
  iteration * ptr = &iters[i-1];
  iteration next = kt.iterate(ptr);
  iters.push_back(next);
  cout << "Done. \n";
  }
};

// ----------------------------------------------------------------------------
// Iterate through the KT equations n times, storing each iteration for comparison
void isobar::print(int n)
{

  cout << endl;
  cout << "Printing the " << std::to_string(n) << "th Iteration..." << endl;
  //Surpress ROOT messages
  gErrorIgnoreLevel = kWarning;

  string name;
  name = "isobar_iteration_" + std::to_string(n);

  // Output to a datfile
  std::ofstream output;
  string namedat = name + ".dat";
  output.open(namedat.c_str());

  vector<double> s;
  vector<double> refx, imfx, ref0, imf0;
  for (int i = 0; i < 60; i++)
  {
    double s_i = sthPi + EPS + double(i) * (.95 - sthPi) / 60.;
    complex<double> fx_i = iters[n].interp_above(s_i);
    complex<double> f0_i = iters[0].interp_above(s_i);

    s.push_back(s_i);
    refx.push_back(real(fx_i)); imfx.push_back(imag(fx_i));
    ref0.push_back(real(f0_i)); imf0.push_back(imag(f0_i));

    output << std::left << setw(15) << s_i << setw(15) << real(fx_i) << setw(15) << imag(fx_i) << endl;
  }
  output.close();

  cout << "Output to: " << namedat << "." << endl;

  //Print the Real part compared to no rescattering
  TCanvas *c = new TCanvas("c", "c");
  TGraph *gRe   = new TGraph(s.size(), &(s[0]), &(refx[0]));
  TGraph *gRe0  = new TGraph(s.size(), &(s[0]), &(ref0[0]));

  string label = std::to_string(n) + " Iterations";
  gRe->SetTitle(label.c_str());
  gRe->SetLineStyle(2);
  gRe->Draw("AL");

  gRe0->Draw("Same");
  gRe0->SetTitle("No Rescattering");

  c->Modified();
  c->BuildLegend();
  string namepdfre = name + "_real.pdf";
  c->Print(namepdfre.c_str());

  delete c, gRe, gRe0;

  //And the Imaginary part
  TCanvas *c2 = new TCanvas("c2", "c2");
  TGraph *gIm   = new TGraph(s.size(), &(s[0]), &(imfx[0]));
  TGraph *gIm0  = new TGraph(s.size(), &(s[0]), &(imf0[0]));

  gIm->SetTitle(label.c_str());
  gIm->SetLineStyle(2);
  gIm->Draw("AL");

  gIm0->Draw("Same");
  gIm0->SetTitle("No Rescattering");

  c2->Modified();
  c2->BuildLegend();
  string namepdfim = name + "_imaginary.pdf";
  c2->Print(namepdfim.c_str());

  cout << "Plot output to: " << namepdfre << ", " << namepdfim << "." << endl;
  delete c2, gIm, gIm0;
};
