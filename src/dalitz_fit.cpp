// Routines to fit a lineshape to another lineshape by minimizing a chi2-type integral over
// the entire dalitz region
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dalitz_fit.hpp"

//-----------------------------------------------------------------------------
// Calculate chi_squared
// Normalized by the area of the Dalitz plot.
double dalitz_fit::chi_squared(const double *par)
{
  dalitz::generate_weights();

  // Update parameters in polynomial expansion at the fitting step
  fit_amp->set_params(n_params, par);

  // Integrate over s
  double s_sum = 0.;
  for (int i = 0; i < dalitz::N_int(); i++)
  {
    double s_i = dalitz::s_abs[i];
    // Integrate over t
    double t_sum = 0.;
    for (int j = 0; j < dalitz::N_int(); j++)
    {
        complex<double> amp_ij, fit_ij;
        double t_ij, ampsqr_ij, fitsqr_ij;

        t_ij = dalitz::t_abs[i][j];

        amp_ij = dalitz::amp->eval(s_i, t_ij);
        ampsqr_ij = abs(amp_ij * amp_ij);

        fit_ij = fit_amp->eval(s_i, t_ij);
        fitsqr_ij = abs(fit_ij * fit_ij);

        double tmp = (ampsqr_ij - fitsqr_ij);

        // Normalize to the center of the dalitz plot
        double normalization = par[0] * abs(dalitz::amp->kinematics.K_jlam(1, 1, dalitz::s_c, dalitz::t_c));

        tmp /= normalization * normalization;

        t_sum += dalitz::t_wgt[i][j] * tmp * tmp;
    };
    s_sum += dalitz::s_wgt[i] * t_sum;
  };

  double chi2 = s_sum / dalitz::dalitz_area();;

  return chi2;
};

// ---------------------------------------------------------------------------
// Minimize the above chi_squared by calling Minuit2 in the ROOT::Math::Minimizer class
void dalitz_fit::extract_params(int eN)
{
  ROOT::Math::Minimizer* minuit = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
  minuit->SetMaxFunctionCalls(10000000);
  minuit->SetTolerance(0.0000001);
  minuit->SetPrintLevel(nError);

  n_params = eN;
  ROOT::Math::Functor fcn(this, &dalitz_fit::chi_squared, n_params);
  minuit->SetFunction(fcn);
  cout << "dalitz_fit: Fitting";

  if (dalitz::amp->kinematics.get_ampName() != "")
  {
    cout << " " + dalitz::amp->kinematics.get_ampName();
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
void dalitz_fit::print_deviation(string filename)
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
    double si = (dalitz::amp->kinematics.smin() + offset) + double(i) * s_step;
    for(int j = 0; j < 100; j++)
    {
     double t_step = (dalitz::amp->kinematics.tmax(si) - 2. * offset - dalitz::amp->kinematics.tmin(si)) / 100.;
     double tij = dalitz::amp->kinematics.tmin(si) + offset  + double(j) * t_step;

     complex<double> ampsqr = dalitz::amp->eval(si, tij);
     complex<double> fitsqr = fit_amp->eval(si, tij);
     double dev = abs(ampsqr / fitsqr) - 1.;
     dev *= 100.;

      output << std::left << setw(15) << dalitz::amp->kinematics.x(si, tij)
                          << setw(15) << dalitz::amp->kinematics.y(si, tij)
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
