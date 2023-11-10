#ifndef INCLUDE_LAMMPS_SYS_HPP
#define INCLUDE_LAMMPS_SYS_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <memory>

#include <btree.hpp>
#include <vec_quat_manipulator.hpp>
#include <boundary_surface.hpp>
#include <atom_array.hpp>
#include <ellipsoid_array.hpp>
#include <bond_array.hpp>
#include <angle_array.hpp>
#include <loop_topology.hpp>

// structure defining bounding box for the system
struct sys_bbox
{
  vec r_min, r_max;
};

// structure containing lengths for spatial system
struct BD_lengths
{
  double r_sphere, r_bdry;
  vec mono_shape, ribo_shape;
};

// structure defining topology switches
struct topo_switches
{
  bool bond;
  bool bending_angle;
  bool twisting_angle;
};


class LAMMPS_sys
{
public:

  // constructor and destructor
  LAMMPS_sys();
  ~LAMMPS_sys();

  // prng seeding
  void prng_seed(int s);

  // set the btree
  void set_btree(btree_state in_state);

  // read the Brownian dynamics lengths
  void read_BD_lengths(std::string lengths_filename);

  // generate a spherical boundary
  void generate_spherical_bdry(double r, double x0, double y0, double z0);

  // read the monomer coordinates and quaternions
  int read_mono_coords(std::string coords_filename, std::string order);
  int read_mono_quats(std::string quats_filename, std::string order);
  // read the ribo coordinates and quaternions
  int read_ribo_coords(std::string coords_filename, std::string order);
  int read_ribo_quats(std::string quats_filename, std::string order);
  // read the bdry coordinates
  int read_bdry_coords(std::string coords_filename, std::string order);

  // write the monomer coordinates and quaternions
  int write_mono_coords(std::string coords_filename, std::string order);
  int write_mono_quats(std::string quats_filename, std::string order);
  // write the ribo coordinates and quaternions
  int write_ribo_coords(std::string coords_filename, std::string order);
  int write_ribo_quats(std::string quats_filename, std::string order);
  // write the bdry coordinates
  int write_bdry_coords(std::string coords_filename, std::string order);

  // apply a mapping to the monomers
  void apply_mono_mapping(std::vector<std::vector<std::array<int,3>>> map);

  // write the monomer coordinates to an xyz for quick visualization
  void write_mono_xyz(std::string data_filename);

  // write the data
  void write_data(std::string data_filename);

  // getters for the total system state
  int get_N_total();
  int get_N_mono();
  int get_N_ribo();
  int get_N_bdry();

  // getters for vectors of properties
  std::vector<vec> get_mono_coords();

  // setters for the subarrays
  void set_N_total(int N);
  void set_N_mono(int N_mono);
  void set_N_ribo(int N_ribo);
  void set_N_bdry(int N_bdry);
  void set_N_total_ellipsoids(int N);
  void set_N_mono_ellipsoids(int N_mono);
  void set_N_ribo_ellipsoids(int N_ribo);

  // get the mono types
  void get_types(int *&t);

  // setters for the total system state
  void set_coords_arr_total(double *&x, std::string order);
  void set_quats_arr_total(double *&q, std::string order);

  // sync the subarrays
  void sync_subarrays();

  // loop topology
  void set_loop_sys_params(loop_sys_params &l_sys_p);
  void initialize_loop_topo(int N_loops);
  void update_loop_topo();
  std::vector<bond> get_loop_bonds();

  // switches for the topology
  void switch_bonds(bool s);
  void switch_bending_angles(bool s);
  void switch_twisting_angles(bool s);
  
private:

  void prepare_test_data();

  // concatenate arrays
  void cat_atom_array(atom_array &in_atoms, atom_array &cat_atoms);
  void cat_ellipsoid_array(ellipsoid_array &in_ellipsoids, ellipsoid_array &cat_ellipsoids);

  // determine the system's topology
  void prepare_topology();
  void set_bonds();
  void set_angles();
  void set_mono_types(int base_type);

  // filters for the topology
  void filter_bonds(int t);
  void filter_angles(int t);

  // prepare the system prior to output
  void finalize_system();
  void merge_system_components();
  void calc_bbox();

  int N_atom_types;
  int N_angle_types, N_bond_types;

  atom_array atoms, mono_atoms, ribo_atoms, bdry_atoms;
  ellipsoid_array ellipsoids, mono_ellipsoids, ribo_ellipsoids;
  bond_array bonds;
  angle_array angles;

  sys_bbox bbox;
  BD_lengths BD_l;
  topo_switches t_s;
  loop_sys_params l_sys_p;
  btree internal_btree;
  boundary_surface b_surf;

  loop_topology loop_topo;

  vec_quat_manipulator vqm;

};

#endif
