#ifndef INCLUDE_ELLIPSOID_ARRAY_HPP
#define INCLUDE_ELLIPSOID_ARRAY_HPP

#include <memory>
#include <iostream>
#include <fstream>
#include <cstring>

#include <vec_quat_manipulator.hpp>

struct ellipsoid
{
  int id;
  vec s;
  quat q;
};

class ellipsoid_array
{
public:

  // constructor and destructor
  ellipsoid_array();
  ~ellipsoid_array();

  // array sizing
  void set_N(int N);
  int get_N();

  // setter and getter for elements
  void set_ellipsoid(int i, ellipsoid e);
  ellipsoid get_ellipsoid(int i);

  // manipulate columns
  void set_min_id(int min_id); // set the minimum atom id and range accordingly [min,min+N)
  void set_shape_all(vec shape); // set the shape of all ellipsoids
  void set_quats(std::vector<quat> qs); // set element-wise orientations
  void set_quats_arr(double *&q, std::string order);
  void set_quat(int i, quat q);
  void normalize_quats();

  // read quats from files
  int read_bin_quats(std::string data_filename, std::string order, bool force_resize);

  // write atom data
  void write(std::fstream &data_file);
  int write_bin(std::string data_filename, std::string order);

private:

  void initialize_ellipsoids();
  void destroy_ellipsoids();

  int N;
  ellipsoid *ellipsoids;
  vec_quat_manipulator vqm;
  
};
#endif
