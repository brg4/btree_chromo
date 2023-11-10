# btree_chromo

## Description

**A C++ program to model theta structures of replicating bacterial chromosomes using binary trees.**

The program is organized in the following manner:

Upon execution a **btree_driver** executes a series of directives stored in a file provided by the user (a full list of directives is provided in the 'Usage' section).

The **btree_driver** contains the following objects.

- **btree** - *a specialized binary tree class representing replicating circular dsDNA in nested theta structures that can only be manipulated using public member functions representing processes physically possible for circular dsDNA*

  	To understand the terminology that will be used to describe replication of circular dsDNA within the program, visualize an analog clock that has been rotated 180 degrees - the *Ori* is at 6 (top) and the *Ter* is at 12 (bottom), monomers are indexed in ascending order (1 -> 12) in the clockwise direction; two replication forks travel in the clockwise (cw, 6->7->8->...) and counter-clockwise (ccw, 6->5->4->...) from the *Ori* towards the *Ter* until they collide. Newly created dsDNA is indexed in ascending order beginning from the end linked to the ccw replication fork to the opposite end of the strand linked to the cw replication fork.

  1) determines topology when system is represented as circular(/theta structure) polymers
  2) counts genome features at nucleotide resolution
  3) prepares coarse-graining into chromosomal loci for contact map calculations
       
- **LAMMPS_sys** - *a class storing the spatial (and other) information of a system of replicating circular dsDNA (monomers of ellipsoidal particles), ribosomes (ellipsoidal particles), and boundary particles (point-like particles) for use in Brownian dynamics simulations using LAMMPS*
       
- **mapper** - *a class that governs the creation of new DNA monomers in the spatial model given an 'initial replication state' and a 'final replication state' using the train-track model of replication*
       
- **LAMMPS_simulator** - *a class that runs Brownian dynamics and accessory routines using LAMMPS to simulate the spatial model*

Member functions of the **btree_driver** control the interactions and exchange of information between these different objects.

By combining a series of directives together a user could perform the following example protocol.

   1) Define a circular chromosome of a specific size.
   2) Read a file containing defined genome features.
   3) Load a model of DnaA replication initiation.
   4) Load the spatial model of a DNA structure in the polymer representation matching the specified size.
   5) Apply transformations stored in a file that represent a series of replication events.
   6) Use the mapper to create new DNA monomers in the spatial model given the change in replication states.
   7) Simulate the spatial model using Brownian dynamics in LAMMPS.
   8) Output the counts of genome features.
   9) Output the coarse-graining for contact map calculations.

   The user now has **A)** an nontrivial replication state, **B)** a matching spatial model, **C)** the copy numbers of genome features in this replication state, and **D)** the coarse-graining necessary to calculate chromosome contact maps for replicating chromosomes.

## Repository directory structure

 - `src` - *source files for program*
 - `include` - *header files for program*
 - `examples` - *examples to demonstrate program*
 - `LAMMPS_src_additions` - *source files to be added to LAMMPS build for DNA polymer model*
 - `LAMMPS_DNA_model` - *files to be included in LAMMPS simulations running DNA polymer model*

## Installation

All dependencies can be installed within a conda virtual environment with the Anaconda python distribution, this is the easiest approach on a local desktop.

   1) Install **OpenMPI** - *I built v4.1.4 from source using GCC-v12.1.0*
   2) Install **fmt Library** - *I built v9.1.0 from source using GCC-v12.1.0*
   3) Follow instructions in '/LAMMPS_src_additions/installation.txt' to make the additions to the LAMMPS source code
   4) Build and install the modified version of **LAMMPS** using cmake with the last command in 'LAMMPS_src_additions/cmake_command.txt' modified to match your build and installation locations (1. cmake command, 2. make, 3. make install) - *I built the 2022/02/17 release using GCC-v12.1.0*
   5) Edit the following variables in the Makefile to match your installations from the previous steps. - *I used environment variables to specify these local installations on my machine in a bash scripts, but you can type them in manually if you would prefer.*
      - **OPENMPI_LDFLAGS** - *OpenMPI library*
      - **FMT_LDFLAGS** - *fmt library*
      - **LAMMPS_LDFLAGS** - *LAMMPS library*
      - **OPENMPI_INCLUDE** - *OpenMPI headers*
      - **FMT_INCLUDE** - *fmt headers*
      - **LAMMPS_INCLUDE** - *LAMMPS headers in source*

   6) `make all` or `make debug`

The executable (**btree_chromo**) will be in `build/apps`.

## Usage

### Getting Started
Prepare a *directives.inp* file containing the directives to be executed by the binary tree program. Lines beginning with '#' are ignored.

Prepare any input files needed for the chosen directives.

Run with: `./btree_chromo (some location)/(some name)_directives.inp`

See `examples` directory for a demonstration:

1) `examples/preparing_chromosome`
2) `examples/querying_chromosome`
3) `examples/preparing_physical_structure`
4) `examples/simulating_chromosome`
5) `examples/simulating_chromosome_with_replication`
6) `examples/simulating_chromosome_with_loops_and_topo`
7) `examples/simulating_ensemble_protocols`

### Program Execution

1) Read the directives
2) Print the directives
3) Parse the directives into commands and parameters
4) Test the validity of the provided parameters
5) Expand any metacommands parsed from metadirectives
6) Use a finite state machine to test the validity of the command sequence
7) Print the final sequence of commands and parameters
8) Execute the command sequence

### Possible Directives (parameters are comma-separated and following ':' when needed)

#### General

 - `terminate` - *terminates the execution immediately*
 - `switch_skip_runs:(T/F)` - *skips all runs/minimizes of the simulator while incrementing timestep for debugging purposes*

#### Controlling Replication State

 - `new_chromo:size` - *initializes an unreplicated chromosome with given size*
 - `input_state:input_file` - *creates state from input_file*
 - `output_state:output_file` - *writes state to output_file*
 - `output_state_at_timestep:output_file` - *execute output_state, but append a modifier to the output_file with the current timestep of the simulator*
 - `print` - *prints binary tree state in terminal*
 - `btree_prng_seed:seed` - *seeds the btree's prng* 
 - `transform:(b)_cw(r_cw)_ccw(r_ccw)` - *applies single transform to branch (b), with replication extents (r\_cw) and (r\_ccw) along clockwise and counter-clockwise directions, respectively*
 - `transforms_file:transforms_file` - *applies transforms stored in transforms_file*
 - `random_transforms:N` - *applies random transforms until (N) units are added or the maximum size is reached*
 
#### Querying Replication State

 - `update_topology` - *solves bond topology of system*
 - `dump_topology:topology_file,idx` - *dumps topology to topology_file with selected indexing convention (idx)*
 - `dump_topology_at_timestep:topology_file,idx` - *execute dump_topology, but append a modifier to the topology_file with the current timestep of the simulator*
 - `dump_fork_partitions:fork_partition_file,idx` - *dumps monomer partitioning about forks to fork_partition_file with selected indexing convention (idx)*
 - `dump_fork_partitions_at_timestep:fork_partition_file,idx` - *execute dump_fork_partitions, but append a modifier to the fork_partition_file with the current timestep of the simulator*
 - `update_CG_map:f_CG` - *update coarse-graining with selected factor (f_CG)*
 - `dump_CG_map:CG_map_file,f_CG,idx` - *dumps CG_map to CG_map_file with selected factor (f_CG) and indexing convention (idx)*
 - `dump_CG_map_at_timestep:CG_map_file,f_CG,idx` - *execute dump_CG_map, but append a modifier to the CG_map_file with the current timestep of the simulator*
 - `regions_file:regions_file,idx` - *reads chromosome regions from regions_file with selected indexing convention (idx)*
 - `dump_regions:regions_count_file,idx` - *updates regions counts given current state and dumps counts to regions_count_file with selected indexing convention (idx)*
 - `dump_regions_at_timestep:regions_count_file` - *execute dump_regions, but append a modifier to the regions_count_file with the current timestep of the simulator*
 
#### Spatial System for Simulations

 - `load_BD_lengths:BD_length_file` - *reads lengths for Brownian dynamics simulation*
 - `spherical_bdry:R,x0,y0,z0` - *generates bdry particles forming a sphere of radius R centered at (x0,y0,z0)*
 - `load_mono_coords:coords_file,order` - *reads binary file with monomer coordinates (doubles) using data ordering convention (row/col)*
 - `load_mono_quats:quats_file,order` - *reads binary file with monomer quaternions (doubles) using data ordering convention (row/col)*
 - `load_ribo_coords:coords_file,order` - *reads binary file with ribosome coordinates (doubles) using data ordering convention (row/col)*
 - `load_ribo_quats:quats_file,order` - *reads binary file with ribosome quaternions (doubles) using data ordering convention (row/col)*
 - `load_bdry_coords:coords_file,order` - *reads binary file with boundary coordinates (doubles) using data ordering convention (row/col)*
 - `write_mono_coords:coords_file,order` - *write binary file with monomer coordinates (doubles) using data ordering convention (row/col)*
 - `write_mono_quats:quats_file,order` - *write binary file with monomer quaternions (doubles) using data ordering convention (row/col)*
 - `write_ribo_coords:coords_file,order` - *write binary file with ribosome coordinates (doubles) using data ordering convention (row/col)*
 - `write_ribo_quats:quats_file,order` - *write binary file with ribosome quaternions (doubles) using data ordering convention (row/col)*
 - `write_bdry_coords:coords_file,order` - *write binary file with boundary coordinates (doubles) using data ordering convention (row/col)*
 - `switch_bonds:(T/F)` - *enable/disable bonds between DNA monomers (default T), must be used prior to 'write_LAMMPS_data_file' to take effect*
 - `switch_bending_angles:(T/F)` - *enable/disable bending angles between DNA monomers (default T), must be used prior to 'write_LAMMPS_data_file' to take effect*
 - `switch_twisting_angles:(T/F)` - *enable/disable twisting angles between DNA monomers (default T), must be used prior to 'write_LAMMPS_data_file' to take effect*
 - `switch_Ori_bdry_attraction:(T/F)` - *enable/disable attraction of Ori bead to bdry using morse/cut potential (default F)*
 - `switch_Ori_pair_repulsion:(T/F)` - *enable/disable repulsion between Ori beads using harmonic/cut potential (default F)*
 - `write_LAMMPS_data:LAMMPS_data_file` - *write a LAMMPS file (data.-) using the current mono, ribo, and bdry coordinates, and the current replication state for the bond/angle topology*
 - `write_mono_xyz:mono_file_xyz` - *write the current monomer coordinates as an .xyz file to load into visualization software*

#### Mapper

 - `set_initial_state` - *set the initial state of the mapper to the current replication state*
 - `set_final_state` - *set the final state of the mapper to the current replication state*
 - `map_replication` - *based on the difference in final and initial replication states, determine new monomer coordinates and add new monomers to the LAMMPS system*
 
#### Simulator

 - `prepare_simulator:log_file` - *initialize MPI and a LAMMPS object that writes its output to log_file, all further commands with 'simulator' in their name will use this LAMMPS object*
 - `simulator_include_file:inc_file` - *executes the 'include' command to run the LAMMPS commands stored in inc_file*
 - `sync_simulator_and_system` - *copies the current simulation state to the LAMMPS_sys object used to control the topology*
 - `clear_simulator` - *execute the 'clear' command to clear the LAMMPS object*
 - `simulator_set_prng_seed:seed` - *seeds the simulator's prng, must be greater than 0*
 - `simulator_set_nProc:nProc` - *sets the number of processors used by the simulator*
 - `simulator_set_DNA_model:DNA_model_dir` - *sets the directory containing the DNA model used by the simulator*
 - `simulator_set_output_details:output_dir,output_label` - *sets the output directory (output_dir) and label (output_label) used for files generated by the simulator*
 - `simulator_set_delta_t:delta_t` - *sets the timestep (delta_t) used for Brownian dynamics within the simulator*
 - `simulator_store_timestep` - *stores the current timestep of the simulator*
 - `simulator_restore_timestep` - *restores the simulator timestep based on the stored timestep*
 - `simulator_increment_timestep:dt` - *increments the timestep by (dt), this must be a non-negative amount*
 - `simulator_reset_timestep:t` - *resets the timestep to (t)*
 - `simulator_reset_prev_dump_timestep:t` - *resets the timestep of the previous dump to (t)*
 - `simulator_read_data:LAMMPS_data_file` - *read a LAMMPS file (data.-) into the simulator*
 - `simulator_minimize_(soft/hard/topoDNA)_(harmonic/FENE):Tfreq` - *run a minimization with the dictated potential while printing thermodynamic information every Tfreq steps*
 - `simulator_run_(soft/hard/topoDNA)_(harmonic/FENE):Nsteps,Tfreq,Dfreq,append_option,skip_option` - *run Brownian dynamics with the dictated potential for Nsteps, while printing thermodynamic information every Tfreq steps and dumping every Dfreq steps - append_option = noappend/append and skip_option = first/skip_first*
 - `simulator_load_loop_params:loop_params_file` - *read a file (loop_params_file) containing the parameters for the looping interactions*
 - `simulator_run_loops:Nloops,Nsteps,Tfreq,Dfreq,append_option,skip_option` - *run Brownian dynamics with the hard/FENE potential for Nsteps with Nloops randomly placed, while printing thermodynamic information every Tfreq steps and dumping every Dfreq steps - append_option = noappend/append and skip_option = first/skip_first*

#### Fused Directives

 - `sys_write_sim_read_LAMMPS_data:LAMMPS_data_file` - *write a LAMMPS file (data.-) with the system, then read the same data file into the simulator*
 - `sys_write_sim_read_LAMMPS_data_at_timestep:LAMMPS_data_file` - *execute sys_write_sim_read_LAMMPS_data, but append a modifier to the LAMMPS_data_file with the current timestep of the simulator*
 - `simulator_relax_progressive:Nsteps,Tfreq` - *run a protocol of 1) minimize_soft_harmonic, 2) run_soft_harmonic, 3) minimize_hard_harmonic, 4) run_hard_harmonic, 5) minimize_soft_FENE to relax the system, where all runs are for Nsteps and thermodynamic information is printed every Tfreq steps*

#### Metadirectives

 - `repeat:N` - *begin a region of directives that will be repeated (N) times*
 - `end_repeat` - *must follow a 'repeat' and terminates the region of directives that will be repeated*
 - `repeat_replicates:min_rep,max_rep,label_padding` - *begin a region of directives that will be repeated for replicates ranging inclusively from (min_rep) to (max_rep), all I/O directives within the region will be modified to include a replicate label of the form '_rep0000X', where (label_padding) specifies the number of zeros*
 - `end_repeat_replicates` - *must follow a 'repeat_replicates' and terminates the region of directives that will be repeated with replicate identifiers*

## Visualization

### Atom Types

 1 - boundary atoms (bdry), 2 - ribosomes (ribo), 3 - DNA monomers (DNA or mono), 4 - Ori monomers (ori), 5 - Ter monomers (ter), 6 - replication fork monomers (fork), 7 - anchor monomers (anchor), 8 - hinge monomers (hinge)

### VMD Instructions

Visualizing the trajectories using VMD requires some extra work to account for the varying atom numbers and atom types.

 1) install LAMMPS plugin for VMD
 2) open TkConsole in VMD
 3) run `set env(LAMMPSDUMMYPOS) {xd,yd,zd}` in the TkConsole, where {xd,yd,zd} is a tuple of the x,y,z coordinates of dummy atoms for systems with varying atom numbers
 4) run `set env(LAMMPSMAXATOMS) Nmax` in the TkConsole, where Nmax is the maximum number of atoms appearing in any frame of the trajectory
 5) run `set env(LAMMPSREMAPFIELDS) {vx=c_id_track,vy=c_type_track}` in the TkConsole, this will remap the fields of "c_id_track" and "c_type_track" to the x and y velocity, respectively, for every frame of the trajectory
 6) Load the trajectory file

The x-velocity (vx) now stores the frame-dependent atom indices and the y-velocity (vy) now stores the frame-dependent atom types for the entire course of the trajectory. Use these fields rather than the default indices and types (which are defined using only the first frame of the trajectory) when making atom selections.

### Ovito Instructions

 Ovito automatically interprets all fields to be frame-dependent.

## Support
brg4@illinois.edu

## Authors and acknowledgment
Benjamin R. Gilbert - brg4@illinois.edu

## Project status
This project is under development.
