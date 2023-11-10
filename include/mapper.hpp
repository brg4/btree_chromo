#ifndef INCLUDE_MAPPER_HPP
#define INCLUDE_MAPPER_HPP

#include <btree.hpp>


class mapper
{
public:

  // constructor and destructor
  mapper();
  ~mapper();

  // setters for states
  void set_initial_state(btree_state st);
  void set_final_state(btree_state st);

  // prepare the mapping
  int prepare_mapping();

  std::vector<std::vector<std::array<int,3>>> get_map();

private:

  // calculate the state difference
  btree_transforms state_diff(btree_state initial_state, btree_state final_state);

  // initialize/destroy the map
  void initialize_map();
  void destroy_map();

  // btrees
  btree_state initial_st, final_st;

  int ***m, *N_new;
  int N_initial, N_final, N_transforms;

};

#endif
