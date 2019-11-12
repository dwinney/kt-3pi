// Auxilary Math Equations
//
// Dependencies: None
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#include "utilities.hpp"

//-----------------------------------------------------------------------------
// Not really math but given int n, outputs right english string
// "1st" vs "2nd" vs "3rd" etc...
// :)
std::string st_nd_rd(int n)
{
  std::string en = std::to_string(n);

  switch (n)
  {
    case 1: en += "st"; break;
    case 2: en += "nd"; break;
    case 3: en += "rd"; break;
    default: en += "th"; break;
  }

  return en;
};

/*******************************************************************************
   Gauss-Legendre integration function, gauleg, from "Numerical Recipes in C"
   (Cambridge Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and
   B.P. Flannery
*******************************************************************************/
#define EPSS 3.0e-11 /* EPSS is the relative precision. */
/******************************************************************************/
void gauleg(double x1, double x2, double x[], double w[], int n)
/*******************************************************************************
   Given the lower and upper limits of integration x1 and x2, and given n, this
   routine returns arrays x[1..n] and w[1..n] of length n, containing the abscissas
   and weights of the Gauss-Legendre n-point quadrature formula.
*******************************************************************************/
{
        int m,j,i;
        double z1,z,xm,xl,pp,p3,p2,p1;
        m=(n+1)/2; /* The roots are symmetric, so we only find half of them. */
        xm=0.5*(x2+x1);
        xl=0.5*(x2-x1);
        for (i=1; i<=m; i++) { /* Loop over the desired roots. */
                z=cos(3.141592654*(i-0.25)/(n+0.5));
                /* Starting with the above approximation to the ith root, we enter */
                /* the main loop of refinement by Newton's method.                 */
                do {
                        p1=1.0;
                        p2=0.0;
                        for (j=1; j<=n; j++) { /* Recurrence to get Legendre polynomial. */
                                p3=p2;
                                p2=p1;
                                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
                        }
                        /* p1 is now the desired Legendre polynomial. We next compute */
                        /* pp, its derivative, by a standard relation involving also  */
                        /* p2, the polynomial of one lower order.                     */
                        pp=n*(z*p1-p2)/(z*z-1.0);
                        z1=z;
                        z=z1-p1/pp; /* Newton's method. */
                } while (fabs(z-z1) > EPSS);
                x[i]=xm-xl*z; /* Scale the root to the desired interval, */
                x[n+1-i]=xm+xl*z; /* and put in its symmetric counterpart.   */
                w[i]=2.0*xl/((1.0-z*z)*pp*pp); /* Compute the weight             */
                w[n+1-i]=w[i]; /* and its symmetric counterpart. */
        }
};

// Utility functions that take in a vector of complex<doubles> and return vector<double>s of the same size
// containing only the real or imaginary parts.
//-----------------------------------------------------------------------------
std::vector<double> vec_real(std::vector<std::complex<double>> fx)
{
  std::vector<double> result;
  for (int i = 0; i < fx.size(); i++)
  {
    result.push_back(real(fx[i]));
  }

  // Quick Error check
  if (result.size() != fx.size())
  {
    std::cout << "vec_real: ERROR Output and Input vector sizes dont match. Quitting... \n";
    std::exit(1);
  }

  return result;
};

std::vector<double> vec_imag(std::vector<std::complex<double>> fx)
{
  std::vector<double> result;
  for (int i = 0; i < fx.size(); i++)
  {
    result.push_back(imag(fx[i]));
  }

  // Quick Error check
  if (result.size() != fx.size())
  {
    std::cout << "vec_imag: ERROR Output and Input vector sizes dont match. Quitting... \n";
    std::exit(1);
  }

  return result;
};

//-----------------------------------------------------------------------------
// Evaluate the interpolation
complex<double> interpolation::operator ()(double s)
{
  complex<double> result(r_inter.Eval(s), i_inter.Eval(s));
  return result;
};

//-----------------------------------------------------------------------------
// Simple function to call ROOT to print a plot
void quick_plot(vector<double> s, vector<complex<double>> fx, string filename)
{
  gErrorIgnoreLevel = kWarning;

  vector<double> refx = vec_real(fx);
  vector<double> imfx = vec_imag(fx);

  TCanvas *c = new TCanvas("c", "c");
  c->Divide(1,2);

  TGraph *gRe   = new TGraph(s.size(), &(s[0]), &(refx[0]));
  TGraph *gIm   = new TGraph(s.size(), &(s[0]), &(imfx[0]));

  TLegend* lRe = new TLegend(0.73,0.77,0.85,0.85);
  TLegend* lIm = new TLegend(0.69,0.77,0.85,0.85);

  c->cd(1);
  gRe->SetTitle(" ");
  gRe->SetLineStyle(2);
  gRe->SetLineWidth(2);
  gRe->SetLineColor(kBlue);
  gRe->Draw("AL");
  lRe->AddEntry(gRe, "REAL PART", "");
  lRe->Draw("same");

  c->cd(2);
  gIm->SetTitle("");
  gIm->SetLineStyle(2);
  gIm->SetLineWidth(2);
  gIm->SetLineColor(kRed);
  gIm->Draw("AL");
  lIm->AddEntry(gIm, "IMAGINARY PART", "");
  lIm->Draw("same");

  c->Modified();

  string namepdf = filename + ".pdf";
  c->Print(namepdf.c_str());
  cout << "Printed to: " << namepdf << std::endl;

  delete c, gRe, gIm;
};

//-----------------------------------------------------------------------------
// Make a pretty PDF file from an input file
// optional parameter is whether to plot using white color scheme to see % difference
void quick_dplot(string file, bool DEV)
{
  gErrorIgnoreLevel = kWarning;

  TCanvas *c = new TCanvas("c", "c");
  TGraph2D *g = new TGraph2D(file.c_str());

  TH2D *h = g->GetHistogram();


  if (DEV == true)
  {
    int NRGBs = 3, NCont = 512;
    gStyle->SetNumberContours(NCont);
    Double_t stops[NRGBs] = { 0.00, 0.50, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 1.00, 1.00 };
    Double_t green[NRGBs] = { 0.00, 1.00, 0.00 };
    Double_t blue[NRGBs]  = { 0.00, 1.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

    h->SetMaximum(3.);
    h->SetMinimum(-3.);
  }
  else
  {
    gStyle->SetPalette(kColorPrintableOnGrey);
  }

  h->SetAxisRange(-1., 1.,"Y");
  h->SetAxisRange(-1., 1.,"X");
  h->Draw("colz");

  c->Modified();
  file.erase(file.end() - 4, file.end());
  file += ".pdf";
  c->Print(file.c_str());

  cout << "Printed to: " << file << std::endl;
  cout << std::endl;

  delete c;
  delete g;
}

//-----------------------------------------------------------------------------
// Print a complex function to a .dat file
void quick_print(vector<double> s, vector<complex<double>> fx, string file)
{
  if (s.size() != fx.size())
  {
    cout << "ERROR: quick_print input vectors not of the same size" << std::endl;
    return;
  }

  string namedat = file + ".dat";
  std::ofstream output;
  output.open(namedat.c_str());

  for (int i = 0; i < s.size(); i++)
  {
    double s_i = s[i];
    complex<double> fx_i = fx[i];

    output << std::left;
    output << std::setw(15) << s_i;
    output << std::setw(15) << real(fx_i);
    output << std::setw(15) << imag(fx_i);
    output << std::setw(15) << abs(fx_i) << std::endl;
  }
  output.close();

  cout << "Output to: " << namedat << std::endl;
};
