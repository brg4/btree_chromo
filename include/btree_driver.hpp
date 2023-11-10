#ifndef INCLUDE_BTREE_DRIVER_HPP
#define INCLUDE_BTREE_DRIVER_HPP

#include <unordered_map>

#include <btree.hpp>
#include <mapper.hpp>
#include <LAMMPS_sys.hpp>
#include <LAMMPS_simulator.hpp>

struct lock
{
  std::string key;
  bool s;
};


class btree_driver //: protected btree
{
public:

  // constructor and destructor
  btree_driver();
  ~btree_driver();

  // read directives from a file
  void read_directives(std::string drctvs_filename);

  // print the directives
  void print_directives();

  // parse the directives
  void parse_directives();

  // validate the parameters in the command sequence
  int validate_command_sequence_parameters();

  // validate the command sequence
  int validate_command_sequence();

  // expand the metacommands
  int expand_metacommands();

  // print the command sequence
  void print_commands();
  
  // execute commands
  int execute_commands();

private:

  // parse a single directive
  void parse_single_directive(std::string drctv, std::string &command, std::vector<std::string> &params);

  // test if metacommands are properly paired
  int test_paired_metacommands(std::string paired_command);
  // compose any metacommands, loops etc.
  void compose_metacommands();
  // expand repeat metacommands
  void expand_repeat_metacommands();
  void expand_repeat_replicates_metacommands();
  void update_replicate_modified_params(std::string &rep_mod, std::string &command, std::vector<std::string> &params);

  // parameter modifiers
  void append_replicate_modifier(std::string &rep_mod, std::string &mod_param);
  void insert_replicate_modifier(std::string &rep_mod, std::string &mod_param);
  std::string get_timestep_modifier();
  void append_timestep_modifier(std::string &ts_mod, std::string &mod_param);
  void insert_timestep_modifier(std::string &ts_mod, std::string &mod_param);
  
  // execute a single command
  int execute_single_command(std::string &command, std::vector<std::string> &params);

  // reset the command locks and updates
  void reset_command_locks_and_updates();
  // set the command requirements
  void prepare_command_requirements();
  lock new_lock(std::string key, bool s);

  // test the command for parameter validity
  int test_command_parameter_validity(std::string &command, std::vector<std::string> &params);
  // test the command for lock state validity
  int test_command_lock_validity(std::string &command);
  // update the lock state following a command
  void update_lock_state_post_command(std::string &command);
  

  //////////////
  // commands // 
  //////////////

  // terminate the directive execution
  int terminate();

  // switch computations on/off
  int switch_skip_runs(std::vector<std::string> &params);

  // create a new chromosome
  int new_chromo(std::vector<std::string> &params);
  
  // input-ouput
  int input_state(std::vector<std::string> &params);
  int output_state(std::vector<std::string> &params);
  int output_state_at_timestep(std::vector<std::string> &params);

  // transforms
  int transforms_file(std::vector<std::string> &params);
  int transform(std::vector<std::string> &params);
  int random_transforms(std::vector<std::string> &params);

  // regions
  int regions_file(std::vector<std::string> &params);
  int dump_regions(std::vector<std::string> &params);
  int dump_regions_at_timestep(std::vector<std::string> &params);

  // topology
  int update_topology();
  int dump_topology(std::vector<std::string> &params);
  int dump_topology_at_timestep(std::vector<std::string> &params);
  int dump_fork_partitions(std::vector<std::string> &params);
  int dump_fork_partitions_at_timestep(std::vector<std::string> &params);

  // coarse-graining
  int update_CG_map(std::vector<std::string> &params);
  int dump_CG_map(std::vector<std::string> &params);
  int dump_CG_map_at_timestep(std::vector<std::string> &params);

  // miscellaneous
  int btree_prng_seed(std::vector<std::string> &params);
  int print_state();

  // LAMMPS system
  // loading coordinates and quaternions
  int load_mono_coords(std::vector<std::string> &params);
  int load_mono_quats(std::vector<std::string> &params);
  int load_ribo_coords(std::vector<std::string> &params);
  int load_ribo_quats(std::vector<std::string> &params);
  int load_bdry_coords(std::vector<std::string> &params);
  // writing coordinates and quaternions
  int write_mono_coords(std::vector<std::string> &params);
  int write_mono_quats(std::vector<std::string> &params);
  int write_ribo_coords(std::vector<std::string> &params);
  int write_ribo_quats(std::vector<std::string> &params);
  int write_bdry_coords(std::vector<std::string> &params);
  // manual boundary specification
  int spherical_bdry(std::vector<std::string> &params);
  // loading BD lengths
  int load_BD_lengths(std::vector<std::string> &params);
  // manipulate system interactions
  int switch_bonds(std::vector<std::string> &params);
  int switch_bending_angles(std::vector<std::string> &params);
  int switch_twisting_angles(std::vector<std::string> &params);
  int switch_extra_potential(std::string extra_pot, std::string s);
  int switch_Ori_bdry_attraction(std::vector<std::string> &params);
  int switch_Ori_pair_repulsion(std::vector<std::string> &params);
  // writing LAMMPS data file
  int write_LAMMPS_data(std::vector<std::string> &params);
  int write_mono_xyz(std::vector<std::string> &params);

  // mapper
  int set_initial_state();
  int set_final_state();
  int map_replication();

  // simulator
  int prepare_simulator(std::vector<std::string> &params);
  int simulator_include_file(std::vector<std::string> &params);
  int sync_simulator_and_system();
  int clear_simulator();
  int simulator_read_data(std::vector<std::string> &params);
  int simulator_set_nProc(std::vector<std::string> &params);
  int simulator_set_prng_seed(std::vector<std::string> &params);
  int simulator_set_DNA_model(std::vector<std::string> &params);
  int simulator_set_output_details(std::vector<std::string> &params);
  int simulator_set_delta_t(std::vector<std::string> &params);
  int simulator_store_timestep();
  int simulator_restore_timestep();
  int simulator_increment_timestep(std::vector<std::string> &params);
  int simulator_reset_prev_dump_timestep(std::vector<std::string> &params);
  int simulator_reset_timestep(std::vector<std::string> &params);
  
  // simulator minimization routines
  template <int SOFT_HARD_TOPO, int HARMONIC_FENE>
  int simulator_minimize(std::vector<std::string> &params);
  
  // simulator run routines
  template <int SOFT_HARD_TOPO, int HARMONIC_FENE>
  int simulator_run(std::vector<std::string> &params);

  // simulator looped DNA routines
  int simulator_load_loop_params(std::vector<std::string> &params);
  int simulator_run_loops(std::vector<std::string> &params);

  // simulator manipulations
  int simulator_expand_bdry_particles(std::vector<std::string> &params);

  // fused commands
  int sys_write_sim_read_LAMMPS_data(std::vector<std::string> &params);
  int sys_write_sim_read_LAMMPS_data_at_timestep(std::vector<std::string> &params);
  int simulator_relax_progressive(std::vector<std::string> &params);

  ///////////////
  // variables //
  ///////////////

  bool skip_runs;

  /////////////
  // objects //
  /////////////

  // internal classes for directive execution
  btree driver_bt;
  mapper driver_mapper;
  LAMMPS_sys driver_lmp_sys;
  LAMMPS_simulator driver_lmp_simulator;

  // internal variables for directive execution
  std::vector<chromo_region> driver_rg; // vector of chromo_regions
  CG_map driver_CG; // coarse-graining map

  // directives, commands, and parameters
  std::vector<std::string> drctvs; // set of directives
  std::vector<std::string> commands; // vector of commands as std::strings
  std::vector<std::vector<std::string>> command_params; // vector command parameters as vectors of std::strings

  // variables to hold requirements for directives
  std::unordered_map<std::string,bool> lock_state; // locks for the command sequence
  std::unordered_map<std::string,std::vector<lock>> lock_tests; // lock requirements for command execution
  std::unordered_map<std::string,std::vector<lock>> lock_updates; // lock updates given successful command execution
  std::unordered_map<std::string,size_t> N_param_reqs; // requirements for the number of parameters

};

#endif
