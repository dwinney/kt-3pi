// General routines for Dalitz plot generation
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dalitz.hpp"

//-----------------------------------------------------------------------------
// Routines related to integrating over dalitz region
void dalitz::set_integration_points(int m)
{
  n = m;
  cout << "Number of integration points changed to " << n << "... \n";
  S_WG_GENERATED = false;
};

// The gauleg function in aux_math.hpp indexes from 1 to N so thats why
// these wrapper functions exist. Also to put them into vectors such that the number of points can change.
void dalitz::generate_s_weights()
{
  double weights[N_int() + 1], abscissas[N_int() + 1];

  double smn = amp->kinematics.smin() + offset;
  double smx = amp->kinematics.smax() - offset;

  gauleg(smn , smx, abscissas, weights, N_int() + 1);

  s_wgt.clear(); s_abs.clear();

  for (int i = 1; i < N_int() + 1; i++)
  {
    s_wgt.push_back(weights[i]);
    s_abs.push_back(abscissas[i]);
  }

  // cout << "param_fit: Gaussian weights for s in the Dalitz Region generated with "
  // << N_int() << " points... \n";
  S_WG_GENERATED = true;
};

// Same but now the bounds in the t variable depend on s.
void dalitz::generate_t_weights(vector<double> s)
{
  // Clear preexisting weights and abcsissas
  t_wgt.clear(); t_abs.clear();
  for (int i = 0; i < N_int(); i++)
  {
    double weights[N_int() + 1], abscissas[N_int() + 1];

    // add 0.001 to push values off the boundary where things may be singular
    double tmn = amp->kinematics.tmin(s[i]) + offset;
    double tmx = amp->kinematics.tmax(s[i]) - offset;

    gauleg(tmn, tmx, abscissas, weights, N_int() + 1);

    vector<double> t_wgt_temp, t_abs_temp;
    for (int j = 1; j < N_int() + 1; j++)
    {
        t_wgt_temp.push_back(weights[j]);
        t_abs_temp.push_back(abscissas[j]);
    }

    t_wgt.push_back(t_wgt_temp);
    t_abs.push_back(t_abs_temp);
  }
  // cout << "param_fit: Gaussian weights for t in the Dalitz Region generated with "
  // << N_int() << " points... \n";
  T_WG_GENERATED = true;
};

void dalitz::generate_weights(){
  if (S_WG_GENERATED == false)
  {
    generate_s_weights();
    generate_t_weights(s_abs);
  }
  else if (T_WG_GENERATED == false)
  {
    generate_t_weights(s_abs);
  }
};

// Calculate the area of the physical Dalitz region
double dalitz::dalitz_area()
{
double t_sum, s_sum, s_i, t_ij;
generate_weights();

s_sum = 0.;
for (int i = 0; i < N_int(); i++)
{
  s_i = s_abs[i];
  t_sum = 0.;
  for (int j = 0; j < N_int(); j++)
  {
    t_ij = t_abs[i][j];
    t_sum +=  t_wgt[i][j] * t_ij;
  }
  s_sum += s_wgt[i] * t_sum;
};
return s_sum;
};

//-----------------------------------------------------------------------------
// Doubly Differential Decay Width
double dalitz::d2Gamma(double s, double t)
{
  complex<double> Fsqr = amp->eval(s, t);
  Fsqr *= Fsqr;
  Fsqr /= 3.;

  return abs(Fsqr / normalization);
};

//-----------------------------------------------------------------------------
// Integrated total decay width
double dalitz::Gamma_total()
{
  generate_weights();

  // Integrate over s
  double s_sum = 0.;
  for (int i = 0; i < N_int(); i++)
  {
    double s_i = s_abs[i];
    // Integrate over t
    double t_sum = 0.;
    for (int j = 0; j < N_int(); j++)
    {
        double t_ij;
        t_ij = t_abs[i][j];
        t_sum += t_wgt[i][j] * d2Gamma(s_i, t_ij);
    };
    s_sum += s_wgt[i] * t_sum;
  };

  return s_sum;
};

//-----------------------------------------------------------------------------
// Make a pretty PDF file from an input file
void dalitz::quick_dalitz(string file)
{
  TCanvas *c = new TCanvas("c", "c");
  TGraph2D *g = new TGraph2D(file.c_str());

  TH2D *h = g->GetHistogram();
  h->SetAxisRange(-1., 1.,"Y");
  h->SetAxisRange(-1., 1.,"X");

  h->Draw("colz");
  gStyle->SetPalette(kColorPrintableOnGrey);

  c->Modified();
  file.erase(file.end() - 4, file.end());
  file += ".pdf";
  c->Print(file.c_str());

  cout << "Printed to: " << file << endl;
  cout << endl;

  delete c;
  delete g;
}

//-----------------------------------------------------------------------------
// Print out a txt file with the Dalitz plot and plot with ROOT
void dalitz::plot(string options)
{
  bool KSF = false, NORMED = false;
  // parse options string
  if (options != "")
  {
  std::istringstream iss(options);
  for (std::string s; iss >> s;)
    {
     if (s == "KSF") {KSF = true;}
     if (s == "normalized") {NORMED = true;}
    }
  }

  // Command Line Message
  cout << "\n";
  cout << "Plotting ";

  if (amp->kinematics.get_ampName() != "")
  {
    cout << "(" << amp->kinematics.get_ampName() << ")";
  }
  else
  {
    if (KSF == true) { cout << "KSF ";}
    cout << "amplitude ";
  }
  if (NORMED == true) { cout << "normalized to center ";}
  cout << "in Dalitz region... \n";

  // File name
  std::string name = amp->kinematics.get_decayParticle();
  if (KSF == true){name += "_KSF";}
  if (NORMED == true){name += "_N";}
  name += "_dalitz_plot.dat";

  std::ofstream output;
  output.open(name.c_str());

  double s_step = (amp->kinematics.smax() - amp->kinematics.smin() - offset)/100.;

  for (int i = 0; i < 100; i++)
  {
    double si = amp->kinematics.smin() + offset + double(i) * s_step;
    for(int j = 0; j < 100; j++)
    {
     double t_step = (amp->kinematics.tmax(si) - offset - amp->kinematics.tmin(si)) / 100.;
     double tij = amp->kinematics.tmin(si) + offset  + double(j) * t_step;

     double gam = d2Gamma(si, tij);

     if (KSF == true){gam /= std::abs(amp->kinematics.Kibble(si, tij));}
     if (NORMED == true){gam /= d2Gamma(s_c, t_c);}
     if (NORMED == true && KSF == true){gam *= std::abs(amp->kinematics.Kibble(s_c, t_c));}

      output << std::left << setw(15) << amp->kinematics.x(si, tij)
                          << setw(15) << amp->kinematics.y(si, tij)
                          << setw(15) << gam << endl;
    }
  }
output.close();

cout << "Output to: " << name << endl;

quick_dalitz(name);
};
