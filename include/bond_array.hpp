#ifndef INCLUDE_BOND_ARRAY_HPP
#define INCLUDE_BOND_ARRAY_HPP

#include <memory>
#include <fstream>

#include <vec_quat_manipulator.hpp>

struct bond
{
  int id, type, i, j;
};

class bond_array
{
public:

  // constructor and destructor
  bond_array();
  ~bond_array();

  // array sizing
  void set_N(int N);
  int get_N();

  // setter for elements
  void set_bond(int i, bond b);

  // getter for elements
  bond get_bond(int i);

  // write bond data
  void write(std::fstream &data_file);

private:

  void initialize_bonds();
  void destroy_bonds();

  int N;
  bond *bonds;
  
};
#endif
