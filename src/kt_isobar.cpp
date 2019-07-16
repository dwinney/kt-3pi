// These objects include the isobars (i.e. the partial waves in each channel)
// they inherite from the amplitude class because they describe the scattering dynamics
// of a single spin and isospin-projected subchannel.
//
// Dependencies: decay_kinematics.hpp, iteration.hpp, kt_equation.hpp
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
  cout << "Storing initial Omnes amplitude... " << endl;

  subtraction_polynomial sub_poly(options.use_conformal);
  vector<subtraction> subtractions;

  for (int n = 0; n < options.max_subs + 1; n++)
  {
    // Need to store values of the Omnes amplitude around the unitarity cut for the
    // analytic continuation of the angular integral.
    vector<double> s;
    vector< complex<double> > above_cut, below_cut;
    complex<double> ab, be;

    for (int i = 0; i < interpolation::N_interp; i++)
    {
      double s_i = sthPi + EPS + double(i) * (omnes::LamOmnes - sthPi) / double(interpolation::N_interp);
      s.push_back(s_i);

      // above the unitarity cut (+ i epsilon)
      ab = omega(s_i, +1) * sub_poly(n, s_i, +1);
      above_cut.push_back(ab);

      // below the cut (- i epsilon)
      be = omega(s_i, -1) * sub_poly(n, s_i, -1);
      below_cut.push_back(be);
    }

  // Make a vector with copies = number of total subtractions
    subtraction not_actually_subtracted(n, s, above_cut, below_cut);
    subtractions.push_back(not_actually_subtracted);
  }

  // If there are previously stored iterations for some reason, clear it when start() is called.
  if (iters.size() != 0)
  {
    iters.clear();
    cout << "isobar: Clearing previously stored Iterations for j = " << spin_proj;
    cout << " and I = " << iso_proj << ". \n";
  }

  iteration zeroth(0, options.max_subs, omega, subtractions);
  iters.push_back(zeroth);

  cout << "Done." << endl;
  cout << endl;
};

// ----------------------------------------------------------------------------
// Iterate through the KT equations n times, storing each iteration for comparison
void isobar::iterate()
{
  start();

  for (int i = 1; i < options.max_iters + 1; i++)
  {
    cout << "Calculating iteration (" << i << "/" << options.max_iters << ")... " << endl;

    // Pass a pointer to the kt_equations
    if (iters.size() >= 1)
    {
      iteration * ptr = &iters[i - 1];
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
// Additionally make a PDF plot comparing the 0th (no rescattering) and the nth iterations.
void isobar::print_iteration(int n, int m)
{

  if (n > iters.size() - 1 ||  m > iters[n].subtractions.size() - 1)
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
  name += "iteration_" + std::to_string(n) + "_" + std::to_string(m);

  // Output to a datfile
  std::ofstream output;
  string namedat = name + ".dat";
  output.open(namedat.c_str());

  vector<double> s;
  vector<double> refx, imfx;
  for (int i = 0; i < 60; i++)
  {
    double s_i = sthPi + EPS + double(i) * (omnes::LamOmnes - sthPi) / 60.;
    complex<double> fx_i =  iters[n].subtractions[m].interp_above(s_i);

    s.push_back(sqrt(s_i));
    refx.push_back(real(fx_i));
    imfx.push_back(imag(fx_i));

    output << std::left << setw(15) << sqrt(s_i) << setw(15) << real(fx_i) << setw(15) << imag(fx_i) << endl;
  }
  output.close();

  cout << "Output to: " << namedat << "." << endl;

  //Print the Real part compared to no rescattering
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
// Set the subtraction coefficients for fitting routines
void isobar::set_params(int n_params, const double *par)
{
  if (n_params != options.max_subs + 1)
  {
    cout << "isobar: Number of parameters and subtrations don't match! Quitting..." << endl;
    exit(1);
  }

  // clear existing coefficients
  coefficients.clear();
  for (int i = 0; i < n_params; i++)
  {
    coefficients.push_back(par[i]);
  }

  // check again for good measure
  if (coefficients.size() != n_params)
  {
  cout << "isobar: Number of parameters and subtrations don't match! Quitting..." << endl;
  exit(1);
  }
};

// ----------------------------------------------------------------------------
// Set the normalization coefficient
void isobar::normalize(double gamma_exp)
{
  cout << "Normalizing isobar amplitude to Gamma_3pi = " << gamma_exp << " MeV..." << endl;

  dalitz<isobar> d_plot(this);
  double gamma = d_plot.Gamma_total();
  double normalization = sqrt(gamma_exp * 1.e-3 / gamma);

  set_params(1, &normalization);
};

void isobar::print_params()
{
    cout << "Printing Subtraction Coefficients... \n";
    cout << "---------------------------------------- \n";
    for (int i = 0; i < options.max_subs + 1; i++)
    {
      cout << std::left <<  setw(20) << "a_" + std::to_string(i) + ":" <<
      setw(20) << coefficients[i] << "\n";
    }
    cout << "---------------------------------------- \n";
    cout << "\n";
};


// ----------------------------------------------------------------------------
// Evaluate the isobar partial wave at some energy s
// Sums subtractions with their coefficients
complex<double> isobar::subtracted_isobar(double s)
{
  if (iters.size() == 0)
  {
      iterate();
  };

  complex<double> result = 0.;
  if (s > sthPi && s <= omega.LamOmnes)
  {
    for (int i = 0; i < options.max_subs + 1; i++)
    {
      result += coefficients[i] * iters.back().subtractions[i].interp_above(s);
    }
  }

  // on the negative real axis (no imaginary part)
  else if (s < sthPi)
  {
    for (int i = 0; i < options.max_subs + 1; i++)
    {
      subtraction_polynomial sub_poly(options.use_conformal);
      result += coefficients[i] * omega(s, 0) * sub_poly(i, s, 0);
    }
  }

  else
  {
    cout << "isobar: Attempting to evaluate outside of interpolation range. Quitting..." << endl;
    exit(1);
  }

  return result;
};

// ----------------------------------------------------------------------------
// Print total isobar including the set coefficients and combining subtractions
void isobar::print()
{
  //Surpress ROOT messages
  gErrorIgnoreLevel = kWarning;

  string name;
  if (kinematics.get_decayParticle() != "")
  {
    name =  kinematics.get_decayParticle() + "_";
  }
  name += "isobar_" + std::to_string(spin_proj) + "_" + std::to_string(iso_proj);

  // Output to a datfile
  std::ofstream output;
  string namedat = name + ".dat";
  output.open(namedat.c_str());

  vector<double> s;
  vector<double> refx, imfx;
  for (int i = 0; i < 60; i++)
  {
    double s_i = (sthPi + EPS) + double(i) * (omnes::LamOmnes - sthPi) / 60.;
    complex<double> fx_i =  subtracted_isobar(s_i);

    s.push_back(sqrt(s_i));
    refx.push_back(real(fx_i));
    imfx.push_back(imag(fx_i));

    output << std::left << setw(15) << sqrt(s_i) << setw(15) << real(fx_i) << setw(15) << imag(fx_i) << endl;
  }
  output.close();

  cout << "Output to: " << namedat << "." << endl;

  //Print the Real part compared to no rescattering
  TCanvas *c = new TCanvas("c", "c");
  c->Divide(1,2);

  TGraph *gRe   = new TGraph(s.size(), &(s[0]), &(refx[0]));
  TGraph *gIm   = new TGraph(s.size(), &(s[0]), &(imfx[0]));

  string label = name;

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
// Evaluate the full isobar amplitude by inserting partial wave in each subchannel
complex<double> isobar::eval(double s, double t)
{
  double u = kinematics.u_man(s,t);
  double zs = kinematics.z_s(s,t), zt = kinematics.z_t(s,t), zu = kinematics.z_u(s,t);

  // complex<double> result = d_hat(1, 1, zs) * subtracted_isobar(s);
  // result += d_hat(1, 1, zt) * subtracted_isobar(t);
  // result += d_hat(1, 1, zu) * subtracted_isobar(u);

  complex<double> result = subtracted_isobar(s);
  result +=  subtracted_isobar(t);
  result +=  subtracted_isobar(u);

  return result;
};
