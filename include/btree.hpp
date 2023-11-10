#ifndef INCLUDE_BTREE_HPP
#define INCLUDE_BTREE_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>
#include <vector>
#include <array>
#include <memory>
#include <random>

struct theta_topo
{
  int start, end, start_link, end_link, mid;
};

struct node
{
  int gen;
  int size;
  int rho_t, rho_cw, rho_ccw;
  bool leaf, complete;
  theta_topo topo;
  node *parent;
  node *left;
  node *right;
};

struct fork_rho
{
  std::string fork;
  int rho_cw, rho_ccw;
};

typedef std::vector<fork_rho> btree_transforms;

struct btree_state
{
  int size;
  //std::vector<fork_rho> fork_rhos;
  btree_transforms transforms;
};

struct chromo_region
{
  std::string name;
  int start, end, count;
};

struct CG_locus
{
  int CG, bCG, start, end;
};

struct CG_leaf
{
  std::string leaf;
  int start, end, mid;
};

struct CG_map
{
  int N_base, N;
  int N_base_CG, N_CG;
  int f_CG;

  std::vector<CG_locus> loci;
  std::vector<CG_leaf> CG_leaves;
};

struct mono_range
{
  bool wrapped;
  int N, ll, ul, mid_ll, mid_ul;
};

struct fork_partition
{
  std::string fork;
  std::vector<mono_range> left_monos, right_monos;
};

class btree
{
public:

  // constructor and destructor
  btree();
  ~btree();

  // reset root
  void reset_root();

  // recursively destroys tree
  void destroy_tree();

  // prng seeding
  void prng_seed(int s);

  // prints entire tree and information
  void print_tree();

  // initialize and branch
  void initialize_tree(int s);
  int branch(std::string loc);

  // prepare and get state
  void prepare_state(btree_state st);
  btree_state get_state();

  // read and write the state
  void write_state(std::string st_filename, btree_state st);
  btree_state read_state(std::string st_filename);
  btree_transforms read_transforms(std::string tr_filename);

  // parse transform of form "(branch)_cw(rho_cw)_ccw(rho_ccw)"
  fork_rho parse_transform(std::string s);

  // apply transformations to btree
  void apply_transforms(btree_transforms tr);
  void single_transform(fork_rho f_r);
  void random_transforms(int r);
  void random_transforms_on_forks(std::vector<std::string> forks, double r);

  // queries about tree state
  int count_total_leaves();
  int count_total_forks();
  int count_completed_forks();
  int count_active_forks();
  int total_size();
  double scaled_size();
  int max_size();
  int single_size();
  int get_generation(std::string loc);

  // solve the theta structure topology
  void solve_topology();
  void dump_topology(std::string topo_filename, int idx);
  theta_topo get_leaf_topo(std::string loc);
  void dump_fork_partitions(std::string fork_partitions_filename, int idx);

  // update the coarse-graining map based on the current state
  CG_map update_CG_map(int f_CG);
  void dump_CG_map(std::string CG_filename, int idx, CG_map &m);

  // read, update, and dump chromo_regions
  std::vector<chromo_region> read_regions(std::string rg_filename, int idx);
  void update_region_counts(std::vector<chromo_region> &c_rs);
  void dump_regions(std::string rg_filename, std::vector<chromo_region> c_rs);

  // get details of tree
  std::vector<std::string> get_completed_forks();
  std::vector<std::string> get_active_forks();
  std::vector<std::string> get_separable_forks();
  std::vector<std::string> get_leaves();

  // create the topologies
  void prepare_bonds(int **&c, int *&t, int &N, int idx);
  void prepare_angles(int **&c, int *&t, int &N, int idx);
  void prepare_types(int *&t, int &N, int base_type);
  
  void foo();

private:

  // recursively destroy tree
  void destroy_tree(node *branch);

  // recursively print branch topology
  void print_branch(node *branch);

  // used for calculating tree state
  int count_leaves(std::string loc);
  int count_total_forks(std::string loc);
  int count_completed_forks(std::string loc);
  int count_active_forks(std::string loc);
  int leaf_counter(node *branch);
  int total_fork_counter(node *branch);
  int completed_fork_counter(node *branch);
  int active_fork_counter(node *branch);
  int branch_size(node *branch);
  int separable_branch_size(node *branch);

  // used for calculating growth
  int grow_at_branch_sym(std::string loc, int r);
  int grow_at_branch_asym(std::string loc, int r_cw, int r_ccw);
  int get_max_growth_cw(node *branch);
  int get_max_growth_ccw(node *branch);
  std::array<int,2> partition_growths_sym(node *branch, int proposed_r_cw, int proposed_r_ccw);

  // get the fork partitions about a fork
  std::vector<fork_partition> get_all_fork_partitions();
  fork_partition get_fork_partition(std::string loc);

  // create centered CG maps per branch
  void centered_CG_map(int &mid, std::vector<CG_locus> &loci, node *branch, int f_CG);

  // traverse leaves and forks to determine identities
  std::vector<std::string> traverse_leaves(std::vector<std::string> leaves, node *branch);
  std::vector<std::string> traverse_completed_forks(std::vector<std::string> forks, int i, node *branch);
  std::vector<std::string> traverse_active_forks(std::vector<std::string> forks, std::string f, node *branch);

  // branch manipulation routines
  node *get_branch(std::string loc);
  void initialize_branch(node *branch, int g, int s);
  void split_branch(node *branch);
  node *parse_dir(node *branch, char d);
  void prepare_Nleaf_state(size_t N_leaves);


  // uniform_real_distribution<> u_dist;
  std::mt19937 rand_eng;
  node *root;
  

};

#endif
