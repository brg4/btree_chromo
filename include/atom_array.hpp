#ifndef INCLUDE_ATOM_ARRAY_HPP
#define INCLUDE_ATOM_ARRAY_HPP

#include <memory>
#include <iostream>
#include <fstream>
#include <cstring>

#include <vec_quat_manipulator.hpp>

struct atom
{
  int id, type, mol_id, ellipsoid_flag;
  double density;
  vec r;
};

class atom_array
{
public:

  // constructor and destructor
  atom_array();
  ~atom_array();

  // array sizing
  void set_N(int N);
  int get_N();

  // manipulate columns
  void reset_ids();
  void set_types(int *&types); // set element-wise types
  void set_type_all(int type); // set type for all atoms
  void set_ellipsoid_flag_all(int ellipsoid_flag); // set ellipsoid flag for all atoms
  void set_mol_ids(int *&mol_ids); // set element-wise mol ids
  void set_mol_id_all(int mol_id); // set single mol id for all atoms
  void set_densities(double *&densities); // set element-wise densities
  void set_density_all(double density); // set single density for all atoms
  void set_coords(std::vector<vec> rs); // set element-wise coordinates
  void set_coords_arr(double *&x, std::string order);
  void set_coord(int i, vec r);

  // read coordinates from files
  int read_bin_coords(std::string data_filename, std::string order, bool force_resize);

  // getter for coordinates
  std::vector<vec> get_coords();

  // setter and getter for elements
  void set_atom(int i, atom a);
  atom get_atom(int i);

  // write atom data
  void write(std::fstream &data_file);
  int write_xyz(std::string data_filename);
  int write_bin(std::string data_filename, std::string order);

private:

  void initialize_atoms();
  void destroy_atoms();

  int N;
  atom *atoms;
  
};
#endif
