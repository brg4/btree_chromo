#ifndef INCLUDE_INTRA_STEP_DISTRIBUTION_HPP
#define INCLUDE_INTRA_STEP_DISTRIBUTION_HPP

#include <iostream>
#include <cmath>

class intra_step_distribution
{
public:

  // constructor and destructor
  intra_step_distribution();
  ~intra_step_distribution();

  // sample a step k
  int get_k(double r, int k_max);

  void print_distribution();

  // prepare a Poisson distribution truncated at k_max and with intensity l
  void poisson_dist(double l, int k_max);
  // prepare a uniform distribution truncated at k_max
  void uniform_dist(int k_max);
  

private:

  // functions for any distribution
  void set_N(int N);
  void reset_w_cw();
  void destroy_w_cw();
  void norm_w_prep_cw();

  // assign uniform weights
  void uniform_w();
  // assign Poisson weights
  void poisson_w(double l);

  
  int N; // total size of distribution
  double *w, *cw; // weights and cumulative weights

};

#endif
