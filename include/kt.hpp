// Self-contained implementation of the GKPY paramterization for low-energy pion scattering.
// Based on: 10.1103/PhysRevD.83.074004
//
// Dependencies: None
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Generic Amplitude Object to store relevant Quantum Numbers for the Decaying Particle
class amplitude
{
protected:
int qn_J, qn_C, qn_P, qn_I, qn_H;
double Mass;

public:
void set_JPC(int j, int p, int c)
{
        qn_J = j; qn_P = p; qn_C = c;
};
void set_Isospin(int i)
{
        qn_I = i;
};
void set_Helicity(int lambda)
{
        qn_H = lambda;
}
void set_Mass(double m)
{
        Mass = m;
};
};
