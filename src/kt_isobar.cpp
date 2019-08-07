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
void isobar::zeroth()
{
  // If there are previously stored iterations for some reason, clear it when start() is called.
  if (subtractions.size() != 0)
  {
    subtractions.clear();
    cout << "isobar: Clearing previously stored subtractions for j = " << spin_proj;
    cout << " and I = " << iso_proj << ". \n";
  }

  subtraction_polynomial sub_poly(options.use_conformal);
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
};

//
// // ----------------------------------------------------------------------------
// // Set the subtraction coefficients for fitting routines
// void isobar::set_params(int n_params, const double *par)
// {
//   if ((n_params - 1) / 2 != options.max_subs)
//   {
//     cout << "isobar: Number of parameters and subtrations don't match! Quitting..." << endl;
//     exit(1);
//   }
//
//   // clear existing coefficients
//   coefficients.clear();
//   normalization = par[0];
//   for (int i = 0; i < (n_params - 1)/2; i++)
//   {
//     complex<double> c_i = par[2*i+1] * exp(xi * par[2*i+2]);
//     coefficients.push_back(c_i);
//   }
//
//   // check again for good measure
//   if (coefficients.size() != (n_params-1) / 2)
//   {
//   cout << "isobar: Number of parameters and subtrations don't match! Quitting..." << endl;
//   exit(1);
//   }
// };
//
// // ----------------------------------------------------------------------------
// // Set the normalization coefficient
// void isobar::normalize(double gamma_exp)
// {
//   cout << endl;
//   cout << "Normalizing total isobar amplitude to Gamma_3pi = " << gamma_exp << " MeV..." << endl;
//
//   dalitz<isobar> d_plot(this);
//   double gamma = d_plot.Gamma_total();
//   normalization = sqrt(gamma_exp * 1.e-3 / gamma);
//
//   cout << "Normalization constant = " << normalization << endl;
//   cout << endl;
// };

// // ----------------------------------------------------------------------------
// // Set second subtraction coefficient to its sum rule values
// void isobar::sum_rule()
// {
//   //check there is two Subtractions
//   if (iters.back().subtractions.size() != 2)
//   {
//     cout << "isobar: need two subtractions to calcualte sum rule values. Quitting... \n";
//     exit(1);
//   }
//
//   cout << "Calculating sum rule value for subtraction constant..." << endl;
//
//   complex<double> b = kt.disp.sum_rule(&iters.back());
//   if (coefficients.size() != 0)
//   {
//     coefficients.clear();
//   };
//   coefficients.push_back(b);
//
//   cout << "Sum rule constant = " << abs(b) <<  " exp(" << arg(b) << "*I)" << endl;
// };

// // ----------------------------------------------------------------------------
// // Print the current stored subtraction coefficients in command line
// void isobar::print_params()
// {
//     cout << "Printing Subtraction Coefficients... \n";
//     cout << "---------------------------------------- \n";
//     cout << std::left << setw(20) << "normalization:" << setw(20) << normalization << endl;
//
//     for (int i = 0; i < options.max_subs; i++)
//     {
//       cout << std::left <<  setw(20) << "|a_" + std::to_string(i) + "|:" <<
//       setw(20) << abs(coefficients[i]) << "\n";
//       cout << std::left <<  setw(20) << "arg a_" + std::to_string(i) + ":" <<
//       setw(20) << arg(coefficients[i]) << "\n";
//     }
//     cout << "---------------------------------------- \n";
//     cout << "\n";
// };
//
//
// // ----------------------------------------------------------------------------
// // Evaluate the isobar partial wave at some energy s
// // Sums subtractions with their coefficients
// complex<double> isobar::subtracted_isobar(double s)
// {
//   if (iters.size() == 0)
//   {
//       cout << "isobar: you need to initialize amplitude and iterate()!. Quitting..." << endl;
//       exit(1);
//   };
//
//   complex<double> result = 0.;
//   if (s > sthPi && s <= omega.LamOmnes)
//   {
//     result += iters.back().subtractions[0].interp_above(s);
//
//     for (int i = 1; i < options.max_subs + 1; i++)
//     {
//       result += coefficients[i-1] * iters.back().subtractions[i].interp_above(s);
//     }
//   }
//
//   else
//   {
//     cout << "isobar: Attempting to evaluate outside of interpolation range. Quitting..." << endl;
//     exit(1);
//   }
//
//   return normalization * result;
// };
//
// // ----------------------------------------------------------------------------
// // Print total isobar including the set coefficients and combining subtractions
// void isobar::print()
// {
//   cout << "Printing isobar with j = "<< std::to_string(spin_proj) << ",";
//   cout << " lambda = " << std::to_string(helicity_proj) << ",";
//   cout << " and I = " << std::to_string(iso_proj);
//   cout << endl;
//
//   //Surpress ROOT messages
//   gErrorIgnoreLevel = kWarning;
//
//   string name;
//   if (kinematics.get_decayParticle() != "")
//   {
//     name =  kinematics.get_decayParticle() + "_";
//   }
//   name += "isobar_" + std::to_string(spin_proj) + "_" + std::to_string(iso_proj);
//
//   // Output to a datfile
//   std::ofstream output;
//   string namedat = name + ".dat";
//   output.open(namedat.c_str());
//
//   vector<double> s;
//   vector<double> refx, imfx;
//   for (int i = 0; i < 60; i++)
//   {
//     double s_i = (sthPi + EPS) + double(i) * (omnes::LamOmnes - sthPi) / 60.;
//     complex<double> fx_i =  subtracted_isobar(s_i);
//
//     s.push_back(sqrt(s_i));
//     refx.push_back(real(fx_i));
//     imfx.push_back(imag(fx_i));
//
//     output << std::left << setw(15) << sqrt(s_i) << setw(15) << real(fx_i) << setw(15) << imag(fx_i);
//     output << setw(15) << abs(fx_i) << endl;
//
//   }
//   output.close();
//
//   cout << "Output to: " << namedat << "." << endl;
//
//   //Print the Real part compared to no rescattering
//   TCanvas *c = new TCanvas("c", "c");
//   c->Divide(1,2);
//
//   TGraph *gRe   = new TGraph(s.size(), &(s[0]), &(refx[0]));
//   TGraph *gIm   = new TGraph(s.size(), &(s[0]), &(imfx[0]));
//
//   string label = name;
//
//   c->cd(1);
//   gRe->SetTitle(label.c_str());
//   gRe->SetLineStyle(2);
//   gRe->SetLineColor(kBlue);
//   gRe->Draw("AL");
//
//   c->cd(2);
//   gIm->SetTitle("Blue = Real part \t \t \t \t \t  Red = Imaginary part");
//   gIm->SetLineStyle(2);
//   gIm->SetLineColor(kRed);
//   gIm->Draw("AL");
//
//   c->Modified();
//   string namepdf = name + ".pdf";
//   c->Print(namepdf.c_str());
//
//   delete c, gRe, gIm;
// };
// // ----------------------------------------------------------------------------
// // Evaluate the full isobar amplitude by inserting partial wave in each subchannel
// complex<double> isobar::eval(double s, double t)
// {
//   double u = kinematics.u_man(s,t);
//   double zs = kinematics.z_s(s,t), zt = kinematics.z_t(s,t), zu = kinematics.z_u(s,t);
//
//   complex<double> result;
//   result  = kinematics.d_hat(spin_proj, helicity_proj, zs) * subtracted_isobar(s);
//   result += kinematics.d_hat(spin_proj, helicity_proj, zt) * subtracted_isobar(t);
//   result += kinematics.d_hat(spin_proj, helicity_proj, zu) * subtracted_isobar(u);
//
//   result *= (2. * spin_proj + 1.);
//   result *= kinematics.K_lambda(helicity_proj, s, t);
//   return result;
// };
