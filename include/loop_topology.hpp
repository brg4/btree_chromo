#ifndef INCLUDE_LOOP_TOPOLOGY_HPP
#define INCLUDE_LOOP_TOPOLOGY_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <btree.hpp>
#include <vec_quat_manipulator.hpp>
#include <loop.hpp>
#include <binding_region.hpp>
#include <intra_step_distribution.hpp>

struct loop_sys_params
{
  int min_dist;
  std::string family;
  double ext_avg;
  int  ext_max;
  double p_unbinding, r_g;
};

struct loop_sim_params
{
  double r_0, k;
  unsigned long freq_loop, freq_topo, dNt_topo;
};

class loop_topology
{
public:

  // constructor and destructor
  loop_topology();
  ~loop_topology();

  // prng seeding
  void prng_seed(int s);

  // set the step distribution
  void set_step_dist(std::string family, double l, int k_max);

  // set the coords
  void set_coords(std::vector<vec> coords);

  // loop topology functions
  void prepare_binding_regions(std::vector<std::string> leaves, std::vector<theta_topo> leaf_topos, int *&t);
  void initialize_loops(int N_loops, int min_dist);
  void update_loops(int ext_max, int min_dist, double p_unbinding, double r_g);
  std::vector<loop> get_loops();
  std::vector<binding_region> get_regions();
  
private:

  // loop topology
  std::vector<loop> loops;
  std::vector<binding_region> regions;
  std::vector<vec> coords;

  // step distribution
  intra_step_distribution step_dist;

  // prng
  std::mt19937 rand_eng;

};

#endif
