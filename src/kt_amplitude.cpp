// Amplitude class that assembles all isobars togethers.
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
    cout << "-> I = 1; j = " << 2*j+1 << "..." << endl;
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

  if (options.max_iters != 0)
  {
    cout << "Starting KT solution..." << endl;
    cout << endl;
  }

  // Iterate up to specified max_iters in the options class
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
// Evaluate the total amplitude by summing the three isobar terms
complex<double> kt_amplitude::eval(double s, double t)
{
  double stu[3] = {s, t, real(kinematics.u_man(s,t))};
  double z_stu[3] = {kinematics.z_s(s,t), kinematics.z_t(s,t), kinematics.z_u(s,t)};

  // Construct amplitude by summing over spins
  // Multiplying by the appropriate prefactors (2j+1) * K_jlam(s,zs) * d_func
  complex<double> result = 0;
  for (int i = 0; 2*i+1 <= options.max_spin; i++)
  {
    isobar * ptr = &iters.back().isobars[i];
    double j = double(ptr->spin_proj);
    double lam = double(ptr->hel_proj);

    // sum each channel
    // TODO: add crossing wigner rotations and sum over lambda^prime
    for (int k = 0; k < 3; k++)
    {
      double x = stu[k], zx = z_stu[k];

      complex<double> temp = (2. * j + 1.);
      temp *= kinematics.barrier_factor(j, lam, x);
      temp *= kinematics.K_jlam(j, lam, x, zx);
      temp *= kinematics.d_hat(j, lam, zx);
      temp *= ptr->subtracted_isobar(x);

      result += temp;
    }
  }

  return normalization * result;
};

// ----------------------------------------------------------------------------
// Set the normalization coefficient by calculating the total decay Width
// and setting it to the input experimental value, gamma_exp
void kt_amplitude::normalize(double gamma_exp)
{
  cout << "Normalizing total amplitude to Gamma_3pi = " << gamma_exp << " MeV..." << endl;

  dalitz d_plot(this);
  double gamma = d_plot.Gamma_total();
  normalization = sqrt(gamma_exp * 1.e-3 / gamma);

  cout << "Normalization constant = " << normalization << endl;
  cout << endl;
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
  vector<complex<double>> fx;
  for (int i = 0; i < 60; i++)
  {
    double s_i = sthPi + EPS + double(i) * (options.interp_cutoff - sthPi) / 60.;
    complex<double> fx_i =  normalization * iters[n].isobars[j].subtractions[m].interp_above(s_i);

    s.push_back(sqrt(s_i));
    fx.push_back(fx_i);

    output << std::left << setw(15) << sqrt(s_i);
    output << setw(15) << real(fx_i) << setw(15) << imag(fx_i);
    output << setw(15) << abs(fx_i) << endl;
  }
  output.close();

  // Plot with ROOT
  quick_plot(s, fx, name);
};

// ----------------------------------------------------------------------------
// Print total isobar including the set coefficients and combining subtractions
void kt_amplitude::print_isobar(int n)
{
  if (n > iters.back().isobars.size() - 1)
  {
    cout << "print_isobar: Cannot print isobar which doesnt exist!" << endl;
    exit(1);
  }

  //Surpress ROOT messages
  gErrorIgnoreLevel = kWarning;

  string name;
  if (kinematics.get_decayParticle() != "")
  {
    name =  kinematics.get_decayParticle() + "_";
  }
  name += "isobar_" + std::to_string(n);

  // Output to a datfile
  std::ofstream output;
  string namedat = name + ".dat";
  output.open(namedat.c_str());

  vector<double> s;
  vector<complex<double>> fx;
  for (int i = 0; i < 60; i++)
  {
    double s_i = (sthPi + EPS) + double(i) * (options.interp_cutoff - sthPi - EPS) / 60.;
    complex<double> fx_i =  normalization * iters.back().isobars[n].subtracted_isobar(s_i);

    s.push_back(sqrt(s_i));
    fx.push_back(fx_i);

    output << std::left << setw(15) << sqrt(s_i) << setw(15) << real(fx_i) << setw(15) << imag(fx_i);
    output << setw(15) << abs(fx_i) << endl;
  }
  output.close();

  cout << "Output to: " << namedat << "." << endl;

  // Plot with ROOT
  quick_plot(s, fx, name);
};
