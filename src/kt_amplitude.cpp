// Amplitude class that assembles all isobars togethers.
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "kt_amplitude.hpp"
// ----------------------------------------------------------------------------
// Display function for all the settings input from the stored kt_options
void kt_amplitude::print_options()
{
  cout << endl;
  cout << "-----------------------------------------------------------" << endl;
  cout << "Using KT equations for ";
  if (kinematics.get_decayParticle() != "")
  {
    cout << kinematics.get_decayParticle() << ":  \n";
  }
  cout << "-> Mass = " << kinematics.get_decayMass() << " GeV, ";
  cout << "J^PC = " << kinematics.print_JPC() << ", \n";
  cout << "-> with max spin j_max = " << options.max_spin << ", \n";
  cout << "-> with ";
  cout << options.max_iters << " iteration(s), \n";
  cout << std::boolalpha << "-> with USE_CONFORMAL = " << options.use_conformal << ", \n";
  cout << "-> with Interpolation up to s = " << options.interp_cutoff << " GeV^2, \n";

  if (options.max_subs == 0)
  {
  cout << "-> and with no free subtraction coefficients." << endl;
  }
  else
  {
    cout << "-> and with " << options.max_subs << " free";
    if (options.use_conformal == true) {cout << " (real)";}
    else {cout << " (complex)";}
    cout << " subtraction coefficient(s):" << endl;
    for (int i = 0; i < options.subIDs.size(); i++)
      {
        cout << "   - " << st_nd_rd(options.subIDs[i][3]) << " order polynomial";
        cout << " on isobar with I = " << options.subIDs[i][0];
        cout << ", j = " << options.subIDs[i][1];
        cout << ", and lambda = " << options.subIDs[i][2] << ". \n";

        if (options.subIDs[i][2] == 0 && kinematics.get_JPC()[1] == -1)
        {
            cout << "CAUTION: Subtraction added for lambda = 0 with P = -1 will not be evaluated." << endl;
        }
        if (options.subIDs[i][1] > options.max_spin)
        {
          cout << "CAUTION: Subtraction added for j > j_max will not be evaluated." << endl;
        }
        if (options.subIDs[i][2] > kinematics.get_totalSpin())
        {
          cout << "CAUTION: Subtraction added for lambda > J will not be evaluated." << endl;
        }
        if ((options.subIDs[i][0] + options.subIDs[i][1]) % 2 != 0)
        {
          cout << "CAUTION: Subtraction added for I + j not even will not be evaluated." << endl;
        }
      }
  }

  cout << "-----------------------------------------------------------" << endl;
};

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

  // need to same bare omnes functions as the zeroth iteration of each isobar
  vector<isobar> bare_omnes;
  for (int j = 0; 2*j+1 <= options.max_spin; j++)
  {
    isobar ith_wave(1, 2*j+1, 1, options, kinematics);

    // isobar::zeroth() stores the bare omnes function with given quantum numbers
    ith_wave.zeroth();
    bare_omnes.push_back(ith_wave);
  }

  // with the vector<isobar> define the zeroth iteration
  iteration zeroth_iter(0, bare_omnes);

  iters.push_back(zeroth_iter);

  cout << "Done." << endl;
};

// ----------------------------------------------------------------------------
// Iterate through the KT equations, storing each iteration for comparison
void kt_amplitude::iterate()
{
  if (options.max_iters != 0)
  {
    cout << endl;
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

      // TODO: Igor's normalization is temp = 1.; NOT (2j+1) change back
      complex<double> temp = xr;
      temp *= kinematics.barrier_factor(j, lam, x);
      temp *= kinematics.K_jlam(j, lam, x, zx);

      temp *= kinematics.d_hat(j, lam, zx);
      temp *= ptr->eval_isobar(x);

      result += temp;
    }
  }

  return normalization * result;
};

// ----------------------------------------------------------------------------
// PARAMETER SETTING

// ----------------------------------------------------------------------------
// Set the normalization coefficient by calculating the total decay Width
// and setting it to the input experimental value, gamma_exp
void kt_amplitude::normalize(double gamma_exp)
{
  cout << "Normalizing total amplitude to Gamma_3pi = " << gamma_exp << " MeV..." << endl;

  dalitz d_plot(this);
  double gamma = d_plot.decay_width();
  normalization = sqrt(gamma_exp * 1.e-3 / gamma);

  cout << "Normalization constant = " << normalization << endl;
  cout << endl;
};

// ----------------------------------------------------------------------------
// Set the subtractions coefficients
void kt_amplitude::set_params(int n, const double * par)
{
  if ( options.use_conformal == false &&
      (n % 2 != 0 || (n / 2) != options.max_subs))
  {
    cout << endl;
    cout << "kt_amplitude::set_params() wrong number of parameters being set! \n";
    cout << "Amplitude has (2 x " << options.max_subs << ") free parameters (modulus and arg)! Quitting... \n";
    exit(0);
  }

  if ( options.use_conformal == true && n != options.max_subs)
  {
    cout << endl;
    cout << "kt_amplitude::set_params() wrong number of parameters being set! \n";
    cout << "Amplitude has " << options.max_subs << " free parameters! Quitting... \n";
    exit(0);
  }

  // loop over the options which specify where and how many subtractions
  int m = 0;
  for (int i = 0; i < options.subIDs.size(); i++)
  {
    // Find isobar belonging to this entry in subIDs
    // i.e. its index in the isobars vector
    int barID;
    if (options.subIDs[i][1] % 2 == 0){barID = options.subIDs[i][1] / 2;}
    else {barID = (options.subIDs[i][1] - 1) / 2;}

    // How many params belong to this isobar
    int ns = options.subIDs[i][3];
    vector<double> iso_par;
    for (int j = 0; j < ns; j++)
    {
      // Save them in pairs if complex or singles if real
      if (options.use_conformal == false)
      {
        iso_par.push_back(par[m + 2*j]);
        iso_par.push_back(par[m + (2*j) + 1]);
      }
      else
      {
        iso_par.push_back(par[m + j]);
      }
    }

    // pass them to the isobar
    iters.back().isobars[barID].set_params(iso_par);

    // continue until all params are allocated
    if (options.use_conformal == false)
    {
      m += 2 * ns;
    }
    else
    {
      m += ns;
    }
  }
};

// ----------------------------------------------------------------------------
// Print the current stored subtraction coefficients in command line
void kt_amplitude::print_params()
{
    cout << endl;
    cout << "-----------------------------------------------------------" << endl;
    cout << "Printing Subtraction Coefficients... \n";

    for (int i = 0; i < options.subIDs.size(); i++)
    {
      int barID;
      if (options.subIDs[i][1] % 2 == 0){barID = options.subIDs[i][1] / 2;}
      else {barID = (options.subIDs[i][1] - 1) / 2;}

      cout << endl;
      cout << "-> j = " << options.subIDs[i][1] << "\n";
      cout << "------------------------------- \n";
      iters.back().isobars[barID].print_params();
    }
    cout << "-----------------------------------------------------------" << endl;
    cout << "\n";
};

// ----------------------------------------------------------------------------
// PLOTTING UTILITIES

// ----------------------------------------------------------------------------
// Print the nth iteration into a dat file.
void kt_amplitude::plot_iteration(int n, int j, int m)
{

  if (n > iters.size() - 1 ||  m > iters[n].isobars[j].subtractions.size() - 1)
  {
      cout << "isobar: Trying to print iteration that doesnt exist. Quitting..." << endl;
      exit(1);
  };

  cout << "Plotting the " + st_nd_rd(n) << " iteration with " << std::to_string(m) << " subtractions..." << endl;

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
void kt_amplitude::plot_isobar(int n)
{
  if (n > iters.back().isobars.size() - 1)
  {
    cout << "print_isobar: Cannot print isobar which doesnt exist!" << endl;
    exit(1);
  }

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
    complex<double> fx_i =  normalization * iters.back().isobars[n].eval_isobar(s_i);

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
