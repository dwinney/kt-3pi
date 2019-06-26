// These are methods for evaluating the KT angular (t) integral (i.e. the angular projection of cross subchannel
// isobars). This is seperated to allow the possibility of testing different methods of evaluating dispersion (s) integral.
//
// Dependencies: kt_iteration.hpp, decay_kinematics.hpp, aux_math.hpp, omnes.hpp
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "kt_ang_integral.hpp"

// ---------------------------------------------------------------------------
// Analytically continued momentum k
complex<double> angular_integral::k(double s)
{
  complex<double> temp1 = sqrt(xr * (s - sthPi) / s);
  temp1 *= sqrt(xr * (b - s)) * sqrt(xr * (a - s));

  return temp1;
};

// ---------------------------------------------------------------------------
// COM Scattering angle in the s-channel scattering subchannel
// Complex in general because of the path of integration needed by analytic continuation
complex<double> angular_integral::z_s(double s, complex<double> t)
{
  complex<double> temp1 = 2. * t + s - mDec * mDec - 3. * mPi * mPi;

  return temp1 / k(s);
};

// ---------------------------------------------------------------------------
// Bounds of integration
complex<double> angular_integral::t_minus(double s)
{
    return (mDec * mDec + 3. * mPi * mPi - s) / 2. - k(s) / 2.;
};

complex<double> angular_integral::t_plus(double s)
{
    return (mDec * mDec + 3. * mPi * mPi - s) / 2. + k(s) / 2.;
};

// ---------------------------------------------------------------------------
// Kinematic kernel inside the angular integral. In general this should depend
// on the helicity and spin of the isobar.
//
// Right now it only has the omega case of Lambda = 1, j = 1, I = I
complex<double> angular_integral::kernel(double s, complex<double> t)
{
  complex<double> zs, temp1;
  zs = z_s(s, t);
  temp1 = (xr - zs * zs);
  temp1 *= pow( xr * (a - s), 1.5);

  return 3. * temp1 / k(s);
};

// ---------------------------------------------------------------------------
// Integration path in the complex plane

// check = 1
// Both limits are real and above the unitarity cut
complex<double> angular_integral::integ_sthPi_a0(double s)
{
    if (s <= sthPi || s >= a0)
    {
      cout << "integ_s0_a0: Integration out of range! Quitting... \n";
      exit(1);
    }

    double w[N_integ + 1], x[N_integ + 1];
    gauleg(real(t_minus(s)), real(t_plus(s)), x, w, N_integ);

    complex<double> sum = 0.;
    for (int i = 1; i < N_integ + 1; i++)
    {
      complex<double> temp = kernel(s, x[i]) * previous->interp_above(x[i]);
      sum += w[i] * temp;
    }

    return sum;
};

// check = 2
// both limits are purely real but one is above the other below
complex<double> angular_integral::integ_a0_a(double s)
{
  if (s < a0|| s > a)
  {
    cout << "integ_a0_a: Integration out of range! Quitting... \n";
    exit(1);
  }

  if (std::abs(t_minus(s) - sthPi + EPS) < 0.0001)
  {
    return 0.;
  }

  double wM[N_integ + 1], xM[N_integ + 1];
  gauleg(real(t_minus(s)), sthPi + EPS, xM, wM, N_integ);

  double wP[N_integ + 1], xP[N_integ + 1];
  gauleg(sthPi + EPS, real(t_plus(s)), xP, wP, N_integ);

  complex<double> sumP = 0., sumM = 0.;
  for(int i = 1; i < N_integ + 1; i++)
  {
    complex<double> tempM = kernel(s, xM[i]) * previous->interp_below(xM[i]);
    complex<double> tempP = kernel(s, xP[i]) * previous->interp_above(xP[i]);

    sumM += wM[i] * tempM;
    sumP += wP[i] * tempP;
  }

  return sumM + sumP;
};

// check = 3
// this is the unphysical region, the bounds of integration are complex to avoid singularities
complex<double> angular_integral::integ_a_b(double s)
{
  if (s < a || s > b)
  {
    cout << "integ_a_b: Integration out of range! Quitting... \n";
    exit(1);
  }

  double w[N_integ + 1], x[N_integ + 1];
  gauleg(0., 1., x, w, N_integ);

  complex<double> sumP = 0., sumM = 0.;
  for(int i = 1; i < N_integ + 1; i++)
  {
    complex<double> z1_i = (1. - x[i]) * t_minus(s) + x[i] * (2. * t_minus(b));
    complex<double> tempM = kernel(s, z1_i) * previous->omega(z1_i, 0);

    tempM *= 2. * t_minus(b) - t_minus(s); // Jacobian
    sumM +=  w[i] * tempM;

    complex<double> z2_i = (1. - x[i]) *  (2. * t_plus(b)) + x[i] * t_plus(s);
    complex<double> tempP = kernel(s, z2_i) * previous->omega(z2_i, 0);

    tempP *= t_plus(s) - 2. * t_plus(b);
    sumP +=  w[i] * tempP;
  }

  return  (sumP + sumM);
};

// check = 4
// t-channel scattering region, limits are real again
complex<double> angular_integral::integ_b(double s)
{
  if (s < b)
  {
    cout << "integ_b: Integration out of range! Quitting... \n";
    exit(1);
  }

  double w[N_integ + 1], x[N_integ + 1];
  gauleg(real(t_minus(s)), real(t_plus(s)), x, w, N_integ);

  complex<double> sum = 0.;
  for (int i = 1; i < N_integ + 1; i++)
  {
    complex<double> temp = kernel(s, x[i]) * previous->omega(x[i], 0);
    sum += w[i] * temp;
  }

  return  sum;
};

// ---------------------------------------------------------------------------
// Angular integral, F_hat
complex<double> angular_integral::operator() (double s)
{
  complex<double> integ;

  if (std::abs(s - sthPi) < 2. * EPS) {
    integ = 2. * previous->interp_above(real(t_minus(sthPi))) * pow(xr * (a - sthPi), 1.5);
  }

  else if (s > sthPi + EPS && s < a0) {
    integ = integ_sthPi_a0(s);
  }

  else if (s >= a0 && s <= a)  {
    integ =  integ_a0_a(s);
  }

  else if (s > a && s < b)  {
    integ = integ_a_b(s);
  }

  else if (s >= b)  {
    integ = integ_b(s);
  }

  else  {
    cout << "AAAAA wtf hella error" << endl;
    exit(1);
  }

  return integ;
};

// ---------------------------------------------------------------------------
// Utitlity function to print out the values of the above integral for testing
void angular_integral::pass_iteration(iteration * prev)
{
    previous = NULL; //reset the pointer for posterity

    previous = prev;
};

// ---------------------------------------------------------------------------
// Utitlity function to print out the values of the above integral for testing
void angular_integral::print(double low, double high)
{
  if (low < sthPi || high > omnes::LamOmnes)
  {
    cout << "angular_integral::print: Error in ploting range. Quitting..." << endl;
    exit(1);
  }

  //Surpress ROOT messages
  gErrorIgnoreLevel = kWarning;

  string name;
  name = "angular_integral";

  // Output to a datfile
  std::ofstream output;
  string namedat = name + ".dat";
  output.open(namedat.c_str());

  vector<double> s;
  vector<double> refx, imfx;

  for (int i = 0; i < 60; i++)
  {
      double s_i = low  + double(i) * (high - low) / 60.;
      complex<double> fx_i = operator()(s_i);

    s.push_back(s_i);
    refx.push_back(real(fx_i)); imfx.push_back(imag(fx_i));

    output << std::left << setw(15) << s_i << setw(15) << real(fx_i) << setw(15) << imag(fx_i) << endl;
  }
  output.close();

  cout << "Output to: " << namedat << "." << endl;

  // Print real part
  TCanvas *c = new TCanvas("c", "c");
  TGraph *gRe   = new TGraph(s.size(), &(s[0]), &(refx[0]));

  gRe->SetTitle("Real Part");
  gRe->Draw("AL");

  c->Modified();
  string namepdfre = name + "_real.pdf";
  c->Print(namepdfre.c_str());

  delete c, gRe;

  //And the Imaginary part
  TCanvas *c2 = new TCanvas("c2", "c2");
  TGraph *gIm   = new TGraph(s.size(), &(s[0]), &(imfx[0]));

  gIm->SetTitle("Imaginary Part");
  gIm->Draw("AL");

  c2->Modified();
  string namepdfim = name + "_imaginary.pdf";
  c2->Print(namepdfim.c_str());

  cout << "Plot output to: " << namepdfre << ", " << namepdfim << "." << endl;
  delete c2, gIm;
};
