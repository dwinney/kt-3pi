#ifndef _OPTS_
#define _OPTS_

#include <vector>
#include <array>

using std::vector;
using std::array;

// Little struct to hold quantities that are global to the KT amplitude
// but not included in decay_kinematics
struct kt_options
{
  // Constructor
  kt_options(){};

  // Copy constructor
  kt_options(const kt_options &old)
  : max_subs(old.max_subs), max_spin(old.max_spin),
    max_iters(old.max_iters),
    use_conformal(old.use_conformal),
    test_angular(old.test_angular),
    interp_cutoff(old.interp_cutoff),
    subIDs(old.subIDs)
  {};

  int max_iters = 0; // Number of total iterations of KT equations
  int max_spin = 0; // Maximal spin projection in sum (not currently implemented)

  bool use_conformal = false; // Whether to use conformal mapping (if TRUE) or standard evaluation (FALSE)
  bool test_angular = false; // Stops the evaluation of the KT equations after calculating the inhomogeneities and printing them to file

  double interp_cutoff = 1.;

  // subtractions
  int max_subs = 0; // total number of free parameters
  vector<array<int, 4>> subIDs; // which isobars are subtrated and how many times

  // User function to add subtractions :)
  void add_subtraction(int iso, int spin, int hel, int n)
  {
    max_subs += n;
    std::array<int,4> sub = {iso, spin, hel, n};
    subIDs.push_back(sub);
  };
};

#endif
