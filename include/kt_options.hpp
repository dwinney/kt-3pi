#ifndef _OPTS_
#define _OPTS_

// Little struct to hold quantities that are global to the KT amplitude
// but not included in decay_kinematics
struct kt_options
{
  int max_subs; // Number of subtractions in dispersion relation
  int max_iters; // Number of total iterations of KT equations
  int max_spin; // Maximal spin projection in sum (not currently implemented)

  bool use_conformal; // Whether to use conformal mapping (if TRUE) or standard evaluation (FALSE)

  kt_options(){};

  kt_options(int max_sub, int max_j, int iter, bool conf)
  : max_subs(max_sub), max_spin(max_j), max_iters(iter), use_conformal(conf)
  {};

  kt_options(const kt_options &old)
  : max_subs(old.max_subs), max_spin(old.max_spin),
    max_iters(old.max_iters),
    use_conformal(old.use_conformal)
  {};
};

#endif
