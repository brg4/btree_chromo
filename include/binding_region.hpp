#ifndef INCLUDE_BINDING_REGION_HPP
#define INCLUDE_BINDING_REGION_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <memory>
#include <random>
#include <unordered_map>

#include <btree.hpp>
#include <vec_quat_manipulator.hpp>

class binding_region
{
public:

  binding_region(std::string leaf, int ll, int ul, int size,
		 bool completed,
		 bool ter_crossing, int mid_ll, int mid_ul);
  ~binding_region();
  
  // get the region size
  int get_size();
  // get the leaf the region belongs to
  std::string get_leaf();
  // get the relative position of a monomer based on the current position, distance, and direction
  int get_relative_monomer_pos(int mono_pos, int dist, int dir);
  // get the relative distance of a monomer
  int get_dist(int mono_pos_i, int mono_pos_j);
  // get the region position from monomer position
  int get_reg_pos(int mono_pos);
  // get the monomer position from the region position
  int get_mono_pos(int reg_pos);

  // randomly select a direction within the region based on the current position
  int select_direction(int mono_pos, int min_dist, double r_d);

  void update_proximities(double r_g, vec a_coord, std::vector<vec> &coords);
  void filter_proximities_near_a(int min_dist, int a_mono_pos);
  std::vector<int> get_and_filter_intra_candidates(int ext_max, int h_mono_pos, int dir);
  std::vector<int> get_inter_candidates();

  // add an unordered_map and an array of indices, then use these in all functions
  void prepare_idx();
  void print_region_map();
  
private:

  // reset the proximities
  void reset_proximities();

  std::string leaf; // leaf that binding region belongs to
  int ll, ul, size; // lower limit and upper limit of indices
  bool completed; // region is a closed circle
  bool ter_crossing; // region contains the crossing at the ter
  int mid_ll, mid_ul; // in positive direction ul->mid_ul->mid_ll->ll

  int *proximities = nullptr; // array containing the proximity results
  
  int *mono_idx = nullptr; // array of monomer indices indexed by region indices
  std::unordered_map<int,int> reg_idx; // map of monomer indices back to region indices
  
  vec_quat_manipulator vqm;
  
};

#endif
