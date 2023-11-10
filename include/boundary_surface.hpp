#ifndef INCLUDE_BOUNDARY_SURFACE
#define INCLUDE_BOUNDARY_SURFACE

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <array>

#include <vec_quat_manipulator.hpp>

struct tri_face
{
  int verts[3];
};

struct edge_map
{
  int vert;
  std::array<int,2> edge;
};

class boundary_surface
{
public:

  // constructor and destructor
  boundary_surface();
  ~boundary_surface();

  // surface initializations
  void unit_icosahedron();
  void unit_tetrahedron();

  // surface preparation
  void generate_sphere(double R, double r);
  
  // surface-wide operations
  void project_to_sphere();
  void scale_coords(double s);
  void interpolate_surface();

  // single face operations
  void interpolate_face(tri_face t_f, std::vector<edge_map> &edge_mapping);

  // getters
  std::vector<vec> get_coords();
  int get_N_verts();

  // write an xyz file for testing
  void write_xyz(std::string data_filename);

private:

  tri_face new_tri_face(int v0, int v1, int v2);

  // edge comparisons
  bool edge_equiv(std::array<int,2> &e0, std::array<int,2> &e1);
  int vert_from_edge(std::array<int,2> &e, std::vector<edge_map> &edge_mapping);
  
  std::vector<tri_face> tri_surf;
  std::vector<vec> coords;

  vec_quat_manipulator vqm;
  
};

#endif
