// Dependencies: None
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "kt_amplitude.hpp"

// ----------------------------------------------------------------------------
// Store initial omnes function for each isobar
// the 0th iteration of KT
void kt_amplitude::start()
{
  // If there are previously stored iterations for some reason, clear it when start() is called.
  if (iters.size() != 0)
  {
    iters.clear();
    cout << "kt_amplitude: Clearing previously stored iterations. \n";
  }

  cout << endl;
  cout << "Storing initial Omnes amplitudes... " << endl;

  //need a vector<isobar> to define the zeroth iteration
  vector<isobar> bare_omnes;
  for (int j = 0; 2*j+1 <= options.max_spin; j++)
  {
    cout << "-> I = 1; j = " << 2*j+1 << "." << endl;
    isobar ith_wave(1, 2*j+1, 1, options, kinematics);

    // isobar::zeroth() stores the bare omnes function with given quantum numbers
    ith_wave.zeroth();
    bare_omnes.push_back(ith_wave);
  }

  iteration zeroth_iter(0, bare_omnes);

  iters.push_back(zeroth_iter);

  cout << "Done." << endl;
  cout << endl;
};

// ----------------------------------------------------------------------------
// Iterate through the KT equations, storing each iteration for comparison
void kt_amplitude::iterate()
{
  start();

  for (int i = 0; i < options.max_iters; i++)
  {
    cout << "Calculating iteration (" << i + 1 << "/" << options.max_iters << ")... " << endl;

    // Pass a pointer of the previous iteration to the kt_equations
    if (iters.size() >= 1)
    {
      iteration * ptr = &iters[i];
      iteration next = kt.iterate(ptr);
      iters.push_back(next); // Save result in the vector iters
    }
    else
    {
      cout << "iterate: Trying to access iteration that is not there. Quitting..." << endl;
      exit(1);
    }

    cout << "Done." << endl;
    cout << endl;
  }
};

// ----------------------------------------------------------------------------
// Print the nth iteration into a dat file.
void kt_amplitude::print_iteration(int n, int j, int m)
{

  if (n > iters.size() - 1 ||  m > iters[n].isobars[j].subtractions.size() - 1)
  {
      cout << "isobar: Trying to print iteration that doesnt exist. Quitting..." << endl;
      exit(1);
  };

  //Surpress ROOT messages
  gErrorIgnoreLevel = kWarning;

  cout << "Printing the " + st_nd_rd(n) << " iteration with " << std::to_string(m) << " subtractions..." << endl;

  string name;
  if (kinematics.get_decayParticle() != "")
  {
    name =  kinematics.get_decayParticle() + "_";
  }
  name += "iteration_" + std::to_string(n) + "_" + std::to_string(j) + "_" + std::to_string(m);

  // Output to a datfile
  std::ofstream output;
  string namedat = name + ".dat";
  output.open(namedat.c_str());

  vector<double> s;
  vector<double> refx, imfx;
  for (int i = 0; i < 60; i++)
  {
    double s_i = sthPi + EPS + double(i) * (1. - sthPi) / 60.;
    complex<double> fx_i =  iters[n].isobars[j].subtractions[m].interp_above(s_i);

    s.push_back(sqrt(s_i));
    refx.push_back(real(fx_i));
    imfx.push_back(imag(fx_i));

    output << std::left << setw(15) << sqrt(s_i) << setw(15) << real(fx_i) << setw(15) << imag(fx_i);
    output << setw(15) << abs(fx_i) << endl;
  }
  output.close();

  cout << "Output to: " << namedat << "." << endl;

  TCanvas *c = new TCanvas("c", "c");
  c->Divide(1,2);

  TGraph *gRe   = new TGraph(s.size(), &(s[0]), &(refx[0]));
  TGraph *gIm   = new TGraph(s.size(), &(s[0]), &(imfx[0]));

  string label = std::to_string(n) + " Iterations " + std::to_string(m) + " Subtractions";

  c->cd(1);
  gRe->SetTitle(label.c_str());
  gRe->SetLineStyle(2);
  gRe->SetLineColor(kBlue);
  gRe->Draw("AL");

  c->cd(2);
  gIm->SetTitle("Blue = Real part \t \t \t \t \t  Red = Imaginary part");
  gIm->SetLineStyle(2);
  gIm->SetLineColor(kRed);
  gIm->Draw("AL");

  c->Modified();
  string namepdf = name + ".pdf";
  c->Print(namepdf.c_str());

  delete c, gRe, gIm;
};

// ----------------------------------------------------------------------------
// Print total isobar including the set coefficients and combining subtractions
void kt_amplitude::print_isobar(int n)
{
  // // cout << "Printing isobar with j = "<< std::to_string(spin_proj) << ",";
  // // cout << " lambda = " << std::to_string(helicity_proj) << ",";
  // // cout << " and I = " << std::to_string(iso_proj);
  // // cout << endl;
  //
  // //Surpress ROOT messages
  // gErrorIgnoreLevel = kWarning;
  //
  // string name;
  // if (kinematics.get_decayParticle() != "")
  // {
  //   name =  kinematics.get_decayParticle() + "_";
  // }
  // // name += "isobar_" + std::to_string(spin_proj) + "_" + std::to_string(iso_proj);
  //
  // // Output to a datfile
  // std::ofstream output;
  // string namedat = name + ".dat";
  // output.open(namedat.c_str());
  //
  // vector<double> s;
  // vector<double> refx, imfx;
  // for (int i = 0; i < 60; i++)
  // {
  //   double s_i = (sthPi + EPS) + double(i) * (omnes::LamOmnes - sthPi) / 60.;
  //   complex<double> fx_i =  subtracted_isobar(s_i);
  //
  //   s.push_back(sqrt(s_i));
  //   refx.push_back(real(fx_i));
  //   imfx.push_back(imag(fx_i));
  //
  //   output << std::left << setw(15) << sqrt(s_i) << setw(15) << real(fx_i) << setw(15) << imag(fx_i);
  //   output << setw(15) << abs(fx_i) << endl;
  //
  // }
  // output.close();
  //
  // cout << "Output to: " << namedat << "." << endl;
  //
  // //Print the Real part compared to no rescattering
  // TCanvas *c = new TCanvas("c", "c");
  // c->Divide(1,2);
  //
  // TGraph *gRe   = new TGraph(s.size(), &(s[0]), &(refx[0]));
  // TGraph *gIm   = new TGraph(s.size(), &(s[0]), &(imfx[0]));
  //
  // string label = name;
  //
  // c->cd(1);
  // gRe->SetTitle(label.c_str());
  // gRe->SetLineStyle(2);
  // gRe->SetLineColor(kBlue);
  // gRe->Draw("AL");
  //
  // c->cd(2);
  // gIm->SetTitle("Blue = Real part \t \t \t \t \t  Red = Imaginary part");
  // gIm->SetLineStyle(2);
  // gIm->SetLineColor(kRed);
  // gIm->Draw("AL");
  //
  // c->Modified();
  // string namepdf = name + ".pdf";
  // c->Print(namepdf.c_str());
  //
  // delete c, gRe, gIm;
};
