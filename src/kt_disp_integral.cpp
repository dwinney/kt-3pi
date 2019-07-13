// These are methods for evaluating the KT dispersion (s) integral.
// This is seperated from the angular integral, and the total KT equations which combine everything
// because of the delicate nature of handleing the singularities in the integrand near pseudo threshold.
//
// Dependencies: kt_iteration.hpp, decay_kinematics.hpp, aux_math.hpp, omnes.hpp
//               angular_integral.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "kt_disp_integral.hpp"

// ---------------------------------------------------------------------------
// The regular (singularity-free) piece of the integrand
// TODO: Here there would be a power of s to the number of subtractions in general
complex<double> dispersion_integral::integrand(int n, double s, int ieps)
{
  complex<double> result = inhom(n, s);
  result *= sin(previous->omega.extrap_phase(s)) / std::abs(previous->omega(s, ieps));

  return result;
};

complex<double> dispersion_integral::integrate(int n, double s, int ieps, double low, double up)
{
  double w[N_integ + 1], x[N_integ + 1];
  gauleg(low, up, x, w, N_integ);

  complex<double> sum = 0.;
  for (int i = 1; i < N_integ + 1; i++)
  {
    complex<double> temp = integrand(n, x[i], ieps);

    // Subtract off F(s) to remove pole at x[i] = s
    temp -= integrand(n, s, ieps);
    temp /=  (x[i] - s - double(ieps) * xi * EPS);

    sum += w[i] * temp;
  }

  return sum / M_PI;
};

complex<double> dispersion_integral::integrate_inf(int n, double s, int ieps, double low)
{
  double w[N_integ + 1], x[N_integ + 1];
  gauleg(0., 1., x, w, N_integ);

  complex<double> sum = 0.;
  for (int i; i < N_integ + 1; i++)
  {
    double sp = low + tan(M_PI * x[i] / 2.);

    complex<double> temp = integrand(n, sp, ieps);
    temp *=  pow(s * xr, options.max_subs) / pow(sp * xr, options.max_subs);
    temp -= integrand(n, s, ieps); // subtract away the pole at s = x[i];
    temp /=  (sp - s - double(ieps) * xi * EPS);
    temp *= s / sp;

    temp *=  (M_PI / 2.) / pow(cos(M_PI * x[i] / 2.), 2.); // jacobian
    sum += w[i] * temp;
  }

  return sum / M_PI;
};

// ---------------------------------------------------------------------------
// Analytical Integral of the subtracted term which regularizes the singularity at s^prime = s
complex<double> dispersion_integral::sp_log(int n, double s, int ieps, double low, double up)
{
  complex<double> result, log_term;
  result = integrand(n, s, ieps);

  if (options.use_conformal == true)
  {
    log_term = log(up - s - double(ieps) * xi * EPS);
    log_term -= log(low - s - double(ieps) * xi * EPS);
  }
  else
  {
    // log_temp = log()
  }

  return result * log_term / M_PI;
};

// ---------------------------------------------------------------------------
// Disperson integral of inhomogeneity
complex<double> dispersion_integral::disperse(int n, double s, int ieps)
{
  complex<double> result;
  // If pseudothreshold is outside of the cutoff integrate without worry
  // no discontinuities
  if (options.use_conformal == true)
  {
    if (LamOmnes < a)
    {
    result = integrate(n, s, ieps, sthPi + EPS, LamOmnes) + sp_log(n, s, ieps, sthPi + EPS, LamOmnes);
    }
    else
    {
      result =  integrate(n, s, ieps, sthPi + EPS, a - interval) + sp_log(n, s, ieps, sthPi + EPS, a - interval);
      result += integrate(n, s, ieps, a + interval, LamOmnes) + sp_log(n, s, ieps, a + interval, LamOmnes);
    }
  }

  else
  {
    result =  integrate(n, s, ieps, sthPi + EPS, a - interval) + sp_log(n, s, ieps, sthPi + EPS, a - interval);
    result += integrate_inf(n, s, ieps, a + interval);
  }

  return result;
};

// ----------------------------------------------------------------------------
// Log term to regularize cut-off point in the conformal mapping scheme
complex<double> dispersion_integral::cutoff_log(int n, double s, int ieps)
{
  complex<double> result = integrand(n, LamOmnes, ieps);
  result /= pow(a * xr - LamOmnes, 1.5);
  result *= log( (LamOmnes - s - xi * double(ieps) * EPS) / (LamOmnes - sthPi));

  return result / M_PI;
};

void dispersion_integral::pass_iteration(iteration * prev)
{
    previous = NULL;

    previous = prev;
    inhom.pass_iteration(prev);
};

// ----------------------------------------------------------------------------
// Evaluate the dispersion integral.
// Method used depends on use_conformal in kt_options
complex<double> dispersion_integral::operator() (int n, double s, int ieps)
{
  // Error check
  if (ieps != -1 && ieps != +1)
  {
    cout << "dispersion_integral: Invalid ieps. Quitting..." << endl;
    exit(1);
  }

  // if testing the angular_integral print file
  if (options.test_angular == true)
  {
    cout << " -> test_angular enabled. Printing inhomogeneity..." << endl;
    angular_test(n);
    cout << "Done." << endl;
    exit(1);
  }

  // evaluate
  if (options.use_conformal == true)
  {
    return disperse(n, s, ieps) - cutoff_log(n, s, ieps);
  }
  else
  {
    return disperse(n, s, ieps);
  }
};

//----------------------------------------------------------------------------
// Print the inhomogeneity
void dispersion_integral::angular_test(int n)
{
  //Surpress ROOT messages
  gErrorIgnoreLevel = kWarning;

  // Output to a datfile
  std::ofstream output;
  string namedat = "inhomogeneity.dat";
  output.open(namedat.c_str());

  vector<double> s;
  vector<double> refx, imfx;
  for (int i = 0; i < 60; i++)
  {
    double s_i = sthPi + 2. * EPS + double(i) * (omnes::LamOmnes - sthPi) / 60.;
    complex<double> fx_i = inhom(n, s_i);

    if (abs(s_i - a) < 0.005)
    {
      continue;
    }

    s.push_back(s_i);
    refx.push_back(real(fx_i));
    imfx.push_back(imag(fx_i));

    output << std::left << setw(15) << s_i << setw(15) << real(fx_i) << setw(15) << imag(fx_i) << endl;
  }
  output.close();

  TCanvas *c = new TCanvas("c", "c");
  c->Divide(1,2);

  TGraph *gRe   = new TGraph(s.size(), &(s[0]), &(refx[0]));
  TGraph *gIm   = new TGraph(s.size(), &(s[0]), &(imfx[0]));

  string label = "Inhomogeneity with " + std::to_string(n) + " Subtractions";

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
  string namepdf = "inhomogeneity.pdf";
  c->Print(namepdf.c_str());

  delete c, gRe, gIm;

  options.test_angular = false;
};
