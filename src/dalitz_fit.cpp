// Routines to fit a lineshape to another lineshape by minimizing a chi2-type integral over
// the entire dalitz region
//
// Dependencies: dalitz.cpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dalitz_fit.hpp"

//-----------------------------------------------------------------------------
// Calculate chi_squared
// Normalized by the area of the Dalitz plot.
template <class T, class F>
double dalitz_fit<T, F>::chi_squared(const double *par)
{
  dalitz<T>::generate_weights();

  // Update parameters in polynomial expansion at the fitting step
  fit_amp->set_params(n_params, par);

  // Integrate over s
  double s_sum = 0.;
  for (int i = 0; i < dalitz<T>::N_int(); i++)
  {
    double s_i = dalitz<T>::s_abs[i];
    // Integrate over t
    double t_sum = 0.;
    for (int j = 0; j < dalitz<T>::N_int(); j++)
    {
        complex<double> amp_ij, fit_ij;
        double t_ij, ampsqr_ij, fitsqr_ij;

        t_ij = dalitz<T>::t_abs[i][j];

        amp_ij = dalitz<T>::amp->eval(s_i, t_ij);
        ampsqr_ij = abs(amp_ij * amp_ij);

        fit_ij = fit_amp->eval(s_i, t_ij);
        fitsqr_ij = abs(fit_ij * fit_ij);

        double tmp = (ampsqr_ij - fitsqr_ij);

        tmp /= abs(dalitz<T>::amp->kinematics.K_lambda(1, dalitz<T>::s_c, dalitz<T>::t_c));
        tmp /= dalitz<T>::amp->error_func(s_i, t_ij);
        tmp /= par[0] * par[0];

        t_sum += dalitz<T>::t_wgt[i][j] * tmp * tmp;
    };
    s_sum += dalitz<T>::s_wgt[i] * t_sum;
  };

  double chi2 = s_sum / dalitz<T>::dalitz_area();;

  return chi2;
};

// ---------------------------------------------------------------------------
// Minimize the above chi_squared by calling Minuit2 in the ROOT::Math::Minimizer class
template <class T, class F>
void dalitz_fit<T, F>::extract_params(int eN)
{
  ROOT::Math::Minimizer* minuit = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
  minuit->SetMaxFunctionCalls(10000000);
  minuit->SetTolerance(0.000001);
  minuit->SetPrintLevel(nError);

  n_params = eN;
  ROOT::Math::Functor fcn(this, &dalitz_fit::chi_squared, n_params);
  minuit->SetFunction(fcn);
  cout << "dalitz_fit: Fitting";

  if (dalitz<T>::amp->kinematics.get_ampName() != "")
  {
    cout << " " + dalitz<T>::amp->kinematics.get_ampName();
  }
   cout << " with " << n_params << " free parameters... \n";

  for (int a = 0; a < n_params ; a++)
  {
    string var_name = "par[" + std::to_string(a) + "]";
    minuit->SetVariable(a, var_name.c_str(), 100., 0.1);
  }

  minuit->Minimize();
  cout << "Minimizing done. " << endl;
  cout << endl;

  double chi2 = minuit->MinValue();
  cout << "sqrt(chi2) = " << sqrt(chi2) * 1.e3 << endl;
  cout << endl;

  fit_amp->set_params(n_params, minuit->X());
};

// ---------------------------------------------------------------------------
// Print the results of a fit as a dalitz plot showing the % deviation from
// my_amp and fit_amp
template <class T, class F>
void dalitz_fit<T, F>::print_deviation(string filename)
{
  gErrorIgnoreLevel = kWarning;

  // Command Line Message
  cout << "Plotting Fit results... \n";

  std::ofstream output;
  if (filename == "")
  {
  filename += "fit_deviation_plot.dat";
  }
  else
  {
    filename += ".dat";
  }
  output.open(filename.c_str());

  double s_step = (smax - smin - 2. * offset) / 100.;

  for (int i = 0; i < 100; i++)
  {
    double si = (dalitz<T>::amp->kinematics.smin() + offset) + double(i) * s_step;
    for(int j = 0; j < 100; j++)
    {
     double t_step = (dalitz<T>::amp->kinematics.tmax(si) - 2. * offset - dalitz<T>::amp->kinematics.tmin(si)) / 100.;
     double tij = dalitz<T>::amp->kinematics.tmin(si) + offset  + double(j) * t_step;

     complex<double> ampsqr = dalitz<T>::amp->eval(si, tij);
     complex<double> fitsqr = fit_amp->eval(si, tij);
     double dev = abs(ampsqr / fitsqr) - 1.;
     dev *= 100.;

      output << std::left << setw(15) << dalitz<T>::amp->kinematics.x(si, tij)
                          << setw(15) << dalitz<T>::amp->kinematics.y(si, tij)
                          << setw(15) << dev << endl;
    }
  }
output.close();

cout << "Output to : " << filename << endl;
cout << endl;

TCanvas *c = new TCanvas("c", "c");
TGraph2D *g = new TGraph2D(filename.c_str());

TH2D *h = g->GetHistogram();

int NRGBs = 3, NCont = 512;
gStyle->SetNumberContours(NCont);
Double_t stops[NRGBs] = { 0.00, 0.50, 1.00 };
Double_t red[NRGBs]   = { 0.00, 1.00, 1.00 };
Double_t green[NRGBs] = { 0.00, 1.00, 0.00 };
Double_t blue[NRGBs]  = { 0.00, 1.00, 0.00 };
TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

h->SetMaximum(3.);
h->SetMinimum(-3.);
h->SetAxisRange(-1., 1.,"Y");
h->SetAxisRange(-1., 1.,"X");
h->SetTitle("");

h->Draw("colz0");
gStyle->SetNumberContours(215);

c->Modified();
filename.erase(filename.end() - 4, filename.end());
filename += ".pdf";
c->Print(filename.c_str());

delete c;
delete g, h;
};
