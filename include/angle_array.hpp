#ifndef INCLUDE_ANGLE_ARRAY_HPP
#define INCLUDE_ANGLE_ARRAY_HPP

#include <memory>
#include <fstream>

#include <vec_quat_manipulator.hpp>

struct angle
{
  int id, type, i, j, k;
};

class angle_array
{
public:

  // constructor and destructor
  angle_array();
  ~angle_array();

  // array sizing
  void set_N(int N);
  int get_N();

  // setter for elements
  void set_angle(int i, angle a);

  // getter for elements
  angle get_angle(int i);

  // write angle data
  void write(std::fstream &data_file);

private:

  void initialize_angles();
  void destroy_angles();

  int N;
  angle *angles;
  
};
#endif
