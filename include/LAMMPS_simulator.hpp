#ifndef INCLUDE_LAMMPS_SIMULATOR_HPP
#define INCLUDE_LAMMPS_SIMULATOR_HPP

// standard library include files
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <memory>
#include <unordered_map>

// OpenMPI include files
#include <mpi.h>

// LAMMPS include files
#include <lammps.h>
#include "domain.h"
#include "group.h"
#include "input.h"
#include "atom.h"
#include "library.h"

// btree_chromo include files
#include <LAMMPS_sys.hpp>

struct thermo_dump_parameters
{
  // bool multi_rep; // multiple replicates requiring replicate number
  bool append, write_first; // append to dump file, write first timestep
  // int rep, rep_padding; // replicate number, padding for replicate label
  int dump_freq, thermo_freq; // dump frequency and thermo frequency
};

class LAMMPS_simulator
{
public:

  // constructor and destructor
  LAMMPS_simulator();
  ~LAMMPS_simulator();

  void LAMMPS_initialize(std::string logfile);
  void LAMMPS_destroy();

  void set_lmp_sys(LAMMPS_sys *lmp_sys);

  // sync the simulator to the system
  void sim_to_sys();
  void sim_to_sys_atom_counts();

  // set nProc for simulator
  void set_nProc(int nProc);
  // set the DNA model
  void set_DNA_model_dir(std::string DNA_model_dir);
  // set output location
  void set_output_details(std::string output_dir, std::string output_file_label);
  // set PRNG seed
  void set_prng_seed(int s);
  // set delta_t
  void set_delta_t(double delta_t);

  // simulation protocol
  void reset_protocol_variables();
  void global_setup();
  void standard_computes();
  void reset_Nt(unsigned long Nt);
  void reset_prev_dump_Nt(unsigned long prev_dump_Nt);
  void command(std::string command);
  void include_file(std::string filename);
  void clear();
  void store_Nt();
  void restore_Nt();

  // read_data
  void read_data(std::string data_file);

  // minimization routines
  void minimize_soft_harmonic(thermo_dump_parameters t_d_p);
  void minimize_hard_harmonic(thermo_dump_parameters t_d_p);
  void minimize_topoDNA_harmonic(thermo_dump_parameters t_d_p);
  void minimize_soft_FENE(thermo_dump_parameters t_d_p);
  void minimize_hard_FENE(thermo_dump_parameters t_d_p);
  void minimize_topoDNA_FENE(thermo_dump_parameters t_d_p);

  // run routines
  void run_soft_harmonic(unsigned long N_steps, thermo_dump_parameters t_d_p);
  void run_hard_harmonic(unsigned long N_steps, thermo_dump_parameters t_d_p);
  void run_topoDNA_harmonic(unsigned long N_steps, thermo_dump_parameters t_d_p);
  void run_soft_FENE(unsigned long N_steps, thermo_dump_parameters t_d_p);
  void run_hard_FENE(unsigned long N_steps, thermo_dump_parameters t_d_p);
  void run_topoDNA_FENE(unsigned long N_steps, thermo_dump_parameters t_d_p);

  // modulate particle size
  void expand_bdry_particles(double ds);
  void reset_bdry_particle_expansion();

  // loop system
  int read_loop_params(std::string loop_param_filename);
  void set_loop_sim_params(loop_sim_params &l_sim_p);
  void update_loop_bonds(bool new_bonds);
  void run_loops(int N_loops, unsigned long N_steps, thermo_dump_parameters t_d_p);

  // switch an extra potential on or off
  void switch_extra_potential(std::string p, bool s);

  unsigned long get_Nt();
  void increment_Nt(unsigned long dNt);
  
private:

  // reset all computes
  void initialize_computes();
  // trigger the state of a compute
  void compute_trigger(std::string compute_label);
  void uncompute(std::string compute_label);

  // reset all dumps
  void initialize_dumps();
  // prepare a .lammpstrj dump
  void prepare_dump(unsigned long N_steps, thermo_dump_parameters &t_d_p);
  void undump(std::string dump_label);

  // setups for runs and minimizes
  void setup_run(unsigned long N_steps, thermo_dump_parameters &t_d_p);
  void setup_minimize(thermo_dump_parameters &t_d_p);
  
  // reset all simulation variables
  void initialize_sim_vars();
  // set a simulation variable to an integer value
  void set_sim_var_int(std::string sim_var, int val);
  void delete_sim_var(std::string sim_var);

  // reset the extra potentials for the simulation
  void initialize_extra_potentials();
  void extra_pots_to_sim_vars();

  // reset the timestep
  void reset_timestep_to_Nt();

  std::unordered_map<std::string,bool> sim_vars; // map storing state of sim_vars
  std::unordered_map<std::string,bool> computes; // map storing state of computes
  std::unordered_map<std::string,bool> dumps; // map storing state of dumps
  std::unordered_map<std::string,bool> extra_pots; // map storing state of extra potentials

  int sim_MPI_initialized, sim_MPI_finalized;
  int sim_MPI_size; // MPI size
  int sim_MPI_rank; // current MPI rank

  unsigned long stored_Nt;
  unsigned long prev_dump_Nt;
  unsigned long Nt;
  int nProc; // number of processors for OpenMP
  int prng_seed; // seed for PRNG within LAMMPS object
  std::string DNA_model_dir, output_dir, output_file_label; // DNA model, output dir, and label for output files
  double delta_t; // timestep size

  // objects

  loop_sim_params l_sim_p;
  
  LAMMPS_NS::LAMMPS *lmp;
  LAMMPS_sys *lmp_sys;
  
  

};

#endif
