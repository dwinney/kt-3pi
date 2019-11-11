// Methods to store and access quantum number information of the decaying particle
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "decay_kinematics.hpp"

//-----------------------------------------------------------------------------
// Set functions

// Set total spin and parity (+-1) and C-parity (+-1)
void decay_kinematics::set_decayJPC(int j, int p, int c)
{
  if (p == 0 || p < -1 || p > 1)
  {
    cout << "Error: in set_decayJPC() P = +1 or -1. Quitting... \n";
    exit(0);
  }

  if (c == 0 || c < -1 || c > 1)
  {
    cout << "Error: in set_decayJPC() C = +1 or -1. Quitting... \n";
    exit(0);
  }

  qn_J = j; qn_P = p; qn_C = c;
};

void decay_kinematics::set_decayIsospin(int i)
{
  if (i > 2 || i < 0)
  {
    cout << "Error: invalid argument in set_decayIsospin(). Quitting... \n";
    exit(0);
  }

  qn_I = i;
};

// s-channel helicity projection of the decay
void decay_kinematics::set_decayHelicity(int lambda)
{
  qn_H = lambda;
};

void decay_kinematics::set_decayMass(double m)
{
  mDec = m;
};

// string to label the name of the decaying particle
void decay_kinematics::set_decayParticle(const char * n)
{
  decay_particle = n;
};

//-----------------------------------------------------------------------------
// Get functions

double decay_kinematics::get_decayMass()
{
  return mDec;
};

int decay_kinematics::get_naturality()
{
  if (qn_P == pow(-1, qn_J))
  {
    return 1;
  }
  else
  {
    return -1;
  }
};

string decay_kinematics::get_decayParticle()
{
  return decay_particle;
};

// return JPC as a three-array 
array<int, 3> decay_kinematics::get_JPC()
{
  array<int,3> result = {qn_J, qn_P, qn_C};
  return result;
};

// Print a nice string to screen with the JPC for example "1-+"
string decay_kinematics::print_JPC()
{
  string result = to_string(qn_J);

  if (qn_P == 1)
  {
    result += "+";
  }
  else if (qn_P == -1)
  {
    result += "-";
  }

  if (qn_C == 1)
  {
    result += "+";
  }
  else if (qn_C == -1)
  {
    result += "-";
  }

  return result;
};

// return just J for conveniece
int decay_kinematics::get_totalSpin()
{
  return qn_J;
}
