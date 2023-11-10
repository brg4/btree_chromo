#include <btree_driver.hpp>

// constructor
btree_driver::btree_driver()
{
  driver_bt.reset_root();
  prepare_command_requirements();
  skip_runs = false;
}


// destructor
btree_driver::~btree_driver()
{
  // std::cout << "btree_driver destructor after destroy" << std::endl;
  driver_bt.destroy_tree();
  driver_rg.clear();
  drctvs.clear();
  // std::cout << "btree_driver destructor after destroy" << std::endl;
} 


// function to read directives from file
void btree_driver::read_directives(std::string drctvs_filename)
{
  std::fstream drctvs_file;
  
  std::string line;

  drctvs_file.open(drctvs_filename, std::ios::in);

  std::cout << "\nREADING DIRECTIVES:\n" << std::endl;
  std::cout << "\t" << drctvs_filename << std::endl;

  if (!drctvs_file.is_open())
    {
      std::cout << "ERROR: file not opened in read_directives" << std::endl;
    }
  else
    {
      while (getline(drctvs_file,line))
	{
	  
	  if ((line.length() > 0) &&
	      (line.at(0) != '#'))
	    {
	      drctvs.push_back(line);
	    }
	      
     	}

      drctvs_file.close();
    }
  
}


// function print all directives
void btree_driver::print_directives()
{
  std::cout << "\n--- PROGRAM DIRECTIVES ---\n" << std::endl;
  for (std::string drctv: drctvs)
    {
      std::cout << drctv << std::endl;
    }
  std::cout << "\n--------------------------\n" << std::endl;
}


// reset the command requirements
void btree_driver::reset_command_locks_and_updates()
{

  // locks
  lock_state["btree_initialized"] = false;
  lock_state["regions_present"] = false;
  lock_state["BD_lengths_present"] = false;
  lock_state["map_initial_present"] = false;
  lock_state["map_final_present"] = false;
  lock_state["simulator_prepared"] = false;
  lock_state["DNA_model"] = false;
  lock_state["output_details"] = false;
  lock_state["delta_t"] = false;
  lock_state["lmp_data_present"] = false;
  lock_state["loop_params_present"] = false;

  // updates
  lock_state["topo_update"] = true;
  lock_state["CG_update"] = true;
  
}


// test the validity of the command given the lock state
int btree_driver::test_command_parameter_validity(std::string &command, std::vector<std::string> &params)
{
  
  // test the number of parameters
  if (N_param_reqs[command] != params.size())
    {
      std::cout << "ERROR: wrong number of parameters for (" << command << ")" << std::endl;
      std::cout << "\t" << params.size() << " were given but " << N_param_reqs[command] << " are required" << std::endl;
      return 1;
    }

  return 0;
}


// validate the parameters in the sequence of commands
int btree_driver::validate_command_sequence_parameters()
{
  std::cout << "\n--- BEGIN COMMAND PARAMETER VALIDATION ---\n" << std::endl;
  
  int e = 0;

  for (size_t i_c=0; i_c<commands.size(); i_c++)
    {
      e += test_command_parameter_validity(commands[i_c],
					   command_params[i_c]);
      if (e > 0) break;
    }

  if (e > 0)
    {
      std::cout << "\terror in command parameters" << std::endl;
    }
  else
    {
      std::cout << "\tvalid command parameters" << std::endl;
    }

  std::cout << "\n--- END COMMAND PARAMETER VALIDATION ---\n" << std::endl;
  return e;
}


// test the validity of the command given the lock state
int btree_driver::test_command_lock_validity(std::string &command)
{

  int e = 0;

  // test the lock state for compatability
  std::string key;
  bool s;
  for (size_t i_lock=0; i_lock<lock_tests[command].size(); i_lock++)
    {
      key = lock_tests[command][i_lock].key;
      s = lock_tests[command][i_lock].s;

      if (lock_state[key] != s)
	{
	  std::cout << "ERROR: incompatible lock state for (" << command << ")" << std::endl;
	  std::cout << "\t" << key << " = " << lock_state[key] << std::endl;
	  e += 1;
	}
    }

  // return the number of errors
  return e;
  
}


// update the locks after the completion of a command
void btree_driver::update_lock_state_post_command(std::string &command)
{
  // update the lock state given the command
  std::string key;
  bool s;
  for (size_t i_lock=0; i_lock<lock_updates[command].size(); i_lock++)
    {
      key = lock_updates[command][i_lock].key;
      s = lock_updates[command][i_lock].s;
      lock_state[key] = s;
    }
}


// create a new lock
lock btree_driver::new_lock(std::string key, bool s)
{
  lock l;
  l.key = key;
  l.s = s;
  return l;
}


// test the validity of the command sequence
int btree_driver::validate_command_sequence()
{

  std::cout << "\n--- BEGIN COMMAND SEQUENCE VALIDATION ---\n" << std::endl;

  reset_command_locks_and_updates();
  
  int e = 0;

  for (size_t i_c=0; i_c<commands.size(); i_c++)
    {
      e += test_command_lock_validity(commands[i_c]);
      update_lock_state_post_command(commands[i_c]);
      if (e > 0) break;
    }

  if (e > 0)
    {
      std::cout << "\terror in command sequence" << std::endl;
    }
  else
    {
      std::cout << "\tvalid command sequence" << std::endl;
    }

  reset_command_locks_and_updates();

  std::cout << "\n--- END COMMAND SEQUENCE VALIDATION ---\n" << std::endl;
  
  return e;
}


// parse all of the directives
void btree_driver::parse_directives()
{
  std::string temp_command;
  std::vector<std::string> temp_params;

  for (std::string drctv : drctvs)
    {
      temp_params.clear();
      parse_single_directive(drctv,temp_command,temp_params);
      commands.push_back(temp_command);
      command_params.push_back(temp_params);
    }
}


// parse a single directive
void btree_driver::parse_single_directive(std::string drctv,
					  std::string &command,
					  std::vector<std::string> &params)
{

  std::string param, temp_params;
  std::string cmd_delim = ":";
  std::string param_delim = ",";
  int delim;

  delim = drctv.find(cmd_delim);

  // get command and parameters
  if (delim != -1)
    {
	  
      command = drctv.substr(0,delim);
      temp_params = drctv.substr(delim+1,drctv.length());

      delim = temp_params.find(param_delim);

      while (delim != -1)
	{
	  param = temp_params.substr(0,delim);
	  params.push_back(param);
	  temp_params.erase(0,delim+1);
	  delim = temp_params.find(param_delim);
	}
      params.push_back(temp_params);
	  
    }
  // get a command without parameters
  else
    {
      command = drctv;
    }
  
}


// print the set of commands
void btree_driver::print_commands()
{
  std::cout << "\n----------------" << std::endl;
  std::cout << "--- COMMANDS ---" << std::endl;
  std::cout << "----------------\n" << std::endl;
  for (size_t i_c=0; i_c<commands.size(); i_c++)
    {
      std::cout << "\nCOMMAND: " << commands[i_c] << std::endl;
      for (size_t i_p=0; i_p<command_params[i_c].size(); i_p++)
	{
	  std::cout << "\tparam_" << i_p
	       << ": " << command_params[i_c][i_p]
	       << std::endl;
	}
    }
  std::cout << "\n--------------\n" << std::endl;
}


// expand the metacommands
int btree_driver::expand_metacommands()
{

  std::cout << "\n--- BEGIN METACOMMAND EXPANSION ---\n" << std::endl;

  int e = 0;

  // test the paired metacommands
  e += test_paired_metacommands("repeat");
  e += test_paired_metacommands("repeat_replicates");

  if (e > 0)
    {
      std::cout << "\terror during metacommand expansion" << std::endl;
    }
  else
    {
      expand_repeat_metacommands();
      expand_repeat_replicates_metacommands();
      std::cout << "\tsuccessful metacommand expansion" << std::endl;
    }

  std::cout << "\n--- END METACOMMAND EXPANSION ---\n" << std::endl;
  
  return e;
}


// test paired metacommands
int btree_driver::test_paired_metacommands(std::string paired_command)
{
  int c = 0;

  for (size_t i_c=0; i_c<commands.size(); i_c++)
    {
      if (commands[i_c] == paired_command)
	{
	  c += 1;
	}
      else if (commands[i_c] == ("end_" + paired_command))
	{
	  c -= 1;
	}
      if (c < 0)
	{
	  break;
	}
    }

  if (c != 0)
    {
      return 1;
    }

  return 0;
}


// expand repeat metacommands
void btree_driver::expand_repeat_metacommands()
{
  int N_repeats;
  size_t i_start=0, i_end=0;
  std::vector<std::string> temp_commands, repeated_commands;
  std::vector<std::vector<std::string>> temp_command_params, repeated_command_params;

  bool repeats_present, start_found, end_found;

  repeats_present = true;

  // loop to find the repeats
  while (repeats_present == true)
    {
      
      // initialize variables for repeat expansion
      start_found = false;
      end_found = false;
      temp_commands.clear();
      temp_command_params.clear();
      repeated_commands.clear();
      repeated_command_params.clear();

      // loop over the commands
      for (size_t i_c=0; i_c<commands.size(); i_c++)
	{
	  if (commands[i_c] == "repeat")
	    {
	      i_start = i_c;
	      start_found = true;
	    }

	  if (commands[i_c] == "end_repeat")
	    {
	      i_end = i_c;
	      end_found = true;
	      break;
	    }
	}
      
      // expand repeats if they are found
      if ((start_found == true) && (end_found == true))
	{

	  // get the number of repeats from the params
	  N_repeats = stoi(command_params[i_start][0]);

	  // add the pre-repeat commands and params
	  for (size_t i_c=0; i_c<i_start; i_c++)
	    {
	      temp_commands.push_back(commands[i_c]);
	      temp_command_params.push_back(command_params[i_c]);
	    }

	  // create the set of repeated commands and params
	  for (size_t i_c=(i_start+1); i_c<i_end; i_c++)
	    {
	      repeated_commands.push_back(commands[i_c]);
	      repeated_command_params.push_back(command_params[i_c]);
	    }

	  // add the repeated commands and params
	  for (int i_rep=0; i_rep<N_repeats; i_rep++)
	    {
	      for (size_t i_c=0; i_c<repeated_commands.size(); i_c++)
		{
		  temp_commands.push_back(repeated_commands[i_c]);
		  temp_command_params.push_back(repeated_command_params[i_c]);
		}
	    }

	  // add the post-repeat commands and params
	  for (size_t i_c=(i_end+1); i_c<commands.size(); i_c++)
	    {
	      temp_commands.push_back(commands[i_c]);
	      temp_command_params.push_back(command_params[i_c]);
	    }

	  // clear the existing commands and params
	  commands.clear();
	  command_params.clear();

	  // copy the temp commands and params
	  commands = temp_commands;
	  command_params = temp_command_params;
	  
	}
      else
	{
	  repeats_present = false;
	}
      
    } // end loop to find repeats
}


// expand repeat metacommands
void btree_driver::expand_repeat_replicates_metacommands()
{
  int min_rep, max_rep, N_reps, padding;
  size_t i_start=0, i_end=0;
  std::vector<std::string> temp_commands, repeated_commands;
  std::vector<std::vector<std::string>> temp_command_params, repeated_command_params;
  std::vector<std::string> replicate_modified_params;

  bool repeats_present, start_found, end_found;

  repeats_present = true;

  // loop to find the repeats
  while (repeats_present == true)
    {
      
      // initialize variables for repeat expansion
      start_found = false;
      end_found = false;
      temp_commands.clear();
      temp_command_params.clear();
      repeated_commands.clear();
      repeated_command_params.clear();

      // loop over the commands
      for (size_t i_c=0; i_c<commands.size(); i_c++)
	{
	  if (commands[i_c] == "repeat_replicates")
	    {
	      i_start = i_c;
	      start_found = true;
	    }

	  if (commands[i_c] == "end_repeat_replicates")
	    {
	      i_end = i_c;
	      end_found = true;
	      break;
	    }
	}
      
      // expand repeats if they are found
      if ((start_found == true) && (end_found == true))
	{

	  // get the minimum, maximum, and padding from the params
	  min_rep = stoi(command_params[i_start][0]);
	  max_rep = stoi(command_params[i_start][1]);
	  padding = stoul(command_params[i_start][2]);

	  // calculate the total number of replicates
	  N_reps = max_rep - min_rep + 1;

	  // add the pre-repeat commands and params
	  for (size_t i_c=0; i_c<i_start; i_c++)
	    {
	      temp_commands.push_back(commands[i_c]);
	      temp_command_params.push_back(command_params[i_c]);
	    }

	  // create the set of repeated commands and params
	  for (size_t i_c=(i_start+1); i_c<i_end; i_c++)
	    {
	      repeated_commands.push_back(commands[i_c]);
	      repeated_command_params.push_back(command_params[i_c]);
	    }

	  // add the repeated commands and params
	  for (int i_rep=0; i_rep<N_reps; i_rep++)
	    {
	      std::string rep_mod = std::to_string(min_rep + i_rep);
	      rep_mod = std::string(padding-std::min<size_t>(padding,rep_mod.length()),'0') + rep_mod;
	      rep_mod = "_rep" + rep_mod;
	      for (size_t i_c=0; i_c<repeated_commands.size(); i_c++)
		{
		  temp_commands.push_back(repeated_commands[i_c]);
		  replicate_modified_params = repeated_command_params[i_c];
		  update_replicate_modified_params(rep_mod,
						   repeated_commands[i_c],
						   replicate_modified_params);
		  temp_command_params.push_back(replicate_modified_params);
		}
	    }

	  // add the post-repeat commands and params
	  for (size_t i_c=(i_end+1); i_c<commands.size(); i_c++)
	    {
	      temp_commands.push_back(commands[i_c]);
	      temp_command_params.push_back(command_params[i_c]);
	    }

	  // clear the existing commands and params
	  commands.clear();
	  command_params.clear();

	  // copy the temp commands and params
	  commands = temp_commands;
	  command_params = temp_command_params;
	  
	}
      else
	{
	  repeats_present = false;
	}
      
    } // end loop to find repeats
}


// update the parameters for commands modified by a replicate number
void btree_driver::update_replicate_modified_params(std::string &rep_mod, std::string &command, std::vector<std::string> &params)
{
  
  if (command == "input_state")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "output_state")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "output_state_at_timestep")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "dump_topology")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "dump_topology_at_timestep")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "dump_fork_partitions")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "dump_fork_partitions_at_timestep")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "dump_CG_map")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "dump_CG_map_at_timestep")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "dump_regions")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "dump_regions_at_timestep")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "load_mono_coords")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "load_mono_quats")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "load_ribo_coords")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "load_ribo_quats")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "load_bdry_coords")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "write_mono_coords")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "write_mono_quats")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "write_ribo_coords")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "write_ribo_quats")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "write_bdry_coords")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "write_LAMMPS_data")
    {
      append_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "write_mono_xyz")
    {
      insert_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "simulator_set_output_details")
    {
      append_replicate_modifier(rep_mod,params[1]);
    }
  else if (command == "simulator_read_data")
    {
      append_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "sys_write_sim_read_LAMMPS_data")
    {
      append_replicate_modifier(rep_mod,params[0]);
    }
  else if (command == "sys_write_sim_read_LAMMPS_data_at_timestep")
    {
      append_replicate_modifier(rep_mod,params[0]);
    }
  
}


// append the replicate modifier
void btree_driver::append_replicate_modifier(std::string &rep_mod, std::string &mod_param)
{
  mod_param = mod_param + rep_mod;
}


// insert the replicate modifier
void btree_driver::insert_replicate_modifier(std::string &rep_mod, std::string &mod_param)
{
  int delim;
  std::string file, file_ext;
  std::string file_ext_delim = ".";

  delim = mod_param.find(file_ext_delim);
  file = mod_param.substr(0,delim);
  file_ext = mod_param.substr(delim,mod_param.length());

  mod_param = file + rep_mod + file_ext; 
}


// get the timestep modifier
std::string btree_driver::get_timestep_modifier()
{
  return "_t" + std::to_string(driver_lmp_simulator.get_Nt());
}

// append the timestep modifier
void btree_driver::append_timestep_modifier(std::string &ts_mod, std::string &mod_param)
{
  mod_param = mod_param + ts_mod;
}


// insert the timestep modifier
void btree_driver::insert_timestep_modifier(std::string &ts_mod, std::string &mod_param)
{
  int delim;
  std::string file, file_ext;
  std::string file_ext_delim = ".";

  delim = mod_param.find(file_ext_delim);
  file = mod_param.substr(0,delim);
  file_ext = mod_param.substr(delim,mod_param.length());

  mod_param = file + ts_mod + file_ext; 
}


// prepare the set of lock tests
void btree_driver::prepare_command_requirements()
{
  lock t_l;
  std::vector<lock> t_ls;


  // terminate
  // number of required parameters
  N_param_reqs["terminate"] = 0;
  // lock tests
  t_ls.clear();
  lock_tests["terminate"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["terminate"] = t_ls;

  // switch_skip_runs
  // number of required parameters
  N_param_reqs["switch_skip_runs"] = 1;
  // lock tests
  t_ls.clear();
  lock_tests["switch_skip_runs"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["switch_skip_runs"] = t_ls;
  

  ////////////////////////////////
  // btree command requirements //
  ////////////////////////////////
  

  // new_chromo
  // number of required parameters
  N_param_reqs["new_chromo"] = 1;
  // lock tests
  t_ls.clear();
  lock_tests["new_chromo"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  t_ls.push_back(new_lock("topo_update",true));
  t_ls.push_back(new_lock("CG_update",true));
  lock_updates["new_chromo"] = t_ls;

  // input_state
  // number of required parameters
  N_param_reqs["input_state"] = 1;
  // lock tests
  t_ls.clear();
  lock_tests["input_state"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  t_ls.push_back(new_lock("topo_update",true));
  t_ls.push_back(new_lock("CG_update",true));
  lock_updates["input_state"] = t_ls;

  // output_state
  // number of required parameters
  N_param_reqs["output_state"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["output_state"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["output_state"] = t_ls;

  // output_state_at_timestep
  // number of required parameters
  N_param_reqs["output_state_at_timestep"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["output_state_at_timestep"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["output_state_at_timestep"] = t_ls;

  // transforms_file
  // number of required parameters
  N_param_reqs["transforms_file"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["transforms_file"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("topo_update",true));
  t_ls.push_back(new_lock("CG_update",true));
  lock_updates["transforms_file"] = t_ls;

  // transform
  // number of required parameters
  N_param_reqs["transform"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["transform"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("topo_update",true));
  t_ls.push_back(new_lock("CG_update",true));
  lock_updates["transform"] = t_ls;

  // random_transforms
  // number of required parameters
  N_param_reqs["random_transforms"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["random_transforms"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("topo_update",true));
  t_ls.push_back(new_lock("CG_update",true));
  lock_updates["random_transforms"] = t_ls;

  // binary_fission
  // number of required parameters
  N_param_reqs["binary_fission"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["binary_fission"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("topo_update",true));
  t_ls.push_back(new_lock("CG_update",true));
  lock_updates["binary_fission"] = t_ls;

  // regions_file
  // number of required parameters
  N_param_reqs["regions_file"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["regions_file"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("regions_present",true));
  lock_updates["regions_file"] = t_ls;

  // dump_regions
  // number of required parameters
  N_param_reqs["dump_regions"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  t_ls.push_back(new_lock("regions_present",true));
  lock_tests["dump_regions"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["dump_regions"] = t_ls;

  // dump_regions_at_timestep
  // number of required parameters
  N_param_reqs["dump_regions_at_timestep"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  t_ls.push_back(new_lock("regions_present",true));
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["dump_regions_at_timestep"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["dump_regions_at_timestep"] = t_ls;

  // update_topology
  // number of required parameters
  N_param_reqs["update_topology"] = 0;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["update_topology"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["update_topology"] = t_ls;

  // dump_topology
  // number of required parameters
  N_param_reqs["dump_topology"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["dump_topology"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["dump_topology"] = t_ls;

  // dump_topology_at_timestep
  // number of required parameters
  N_param_reqs["dump_topology_at_timestep"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["dump_topology_at_timestep"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["dump_topology_at_timestep"] = t_ls;

  // dump_fork_partitions
  // number of required parameters
  N_param_reqs["dump_fork_partitions"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["dump_fork_partitions"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["dump_fork_partitions"] = t_ls;

  // dump_fork_partitions_at_timestep
  // number of required parameters
  N_param_reqs["dump_fork_partitions_at_timestep"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["dump_fork_partitions_at_timestep"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["dump_fork_partitions_at_timestep"] = t_ls;

  // update_CG_map
  // number of required parameters
  N_param_reqs["update_CG_map"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["update_CG_map"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["update_CG_map"] = t_ls;

  // dump_CG_map
  // number of required parameters
  N_param_reqs["dump_CG_map"] = 3;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["dump_CG_map"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["dump_CG_map"] = t_ls;

  // dump_CG_map_at_timestep
  // number of required parameters
  N_param_reqs["dump_CG_map_at_timestep"] = 3;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["dump_CG_map_at_timestep"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["dump_CG_map_at_timestep"] = t_ls;

  // btree_prng_seed
  // number of required parameters
  N_param_reqs["btree_prng_seed"] = 1;
  // lock tests
  t_ls.clear();
  lock_tests["btree_prng_seed"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["btree_prng_seed"] = t_ls;

  // print
  // number of required parameters
  N_param_reqs["print"] = 0;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["print"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["print"] = t_ls;

  
  /////////////////////////////////////
  // LAMMPS_sys command requirements //
  /////////////////////////////////////

  
  // load_BD_lengths
  // number of required parameters
  N_param_reqs["load_BD_lengths"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["load_BD_lengths"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("BD_lengths_present",true));
  lock_updates["load_BD_lengths"] = t_ls;

  // load_mono_coords
  // number of required parameters
  N_param_reqs["load_mono_coords"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["load_mono_coords"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["load_mono_coords"] = t_ls;

  // load_mono_quats
  // number of required parameters
  N_param_reqs["load_mono_quats"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["load_mono_quats"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["load_mono_quats"] = t_ls;

  // load_ribo_coords
  // number of required parameters
  N_param_reqs["load_ribo_coords"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["load_ribo_coords"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["load_ribo_coords"] = t_ls;

  // load_ribo_quats
  // number of required parameters
  N_param_reqs["load_ribo_quats"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["load_ribo_quats"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["load_ribo_quats"] = t_ls;

  // load_bdry_coords
  // number of required parameters
  N_param_reqs["load_bdry_coords"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["load_bdry_coords"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["load_bdry_coords"] = t_ls;

  // write_mono_coords
  // number of required parameters
  N_param_reqs["write_mono_coords"] = 2;
  // lock tests
  t_ls.clear();
  lock_tests["write_mono_coords"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["write_mono_coords"] = t_ls;

  // write_mono_quats
  // number of required parameters
  N_param_reqs["write_mono_quats"] = 2;
  // lock tests
  t_ls.clear();
  lock_tests["write_mono_quats"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["write_mono_quats"] = t_ls;

  // write_ribo_coords
  // number of required parameters
  N_param_reqs["write_ribo_coords"] = 2;
  // lock tests
  t_ls.clear();
  lock_tests["write_ribo_coords"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["write_ribo_coords"] = t_ls;

  // write_ribo_quats
  // number of required parameters
  N_param_reqs["write_ribo_quats"] = 2;
  // lock tests
  t_ls.clear();
  lock_tests["write_ribo_quats"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["write_ribo_quats"] = t_ls;

  // write_bdry_coords
  // number of required parameters
  N_param_reqs["write_bdry_coords"] = 2;
  // lock tests
  t_ls.clear();
  lock_tests["write_bdry_coords"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["write_bdry_coords"] = t_ls;

  // write_LAMMPS_data
  // number of required parameters
  N_param_reqs["write_LAMMPS_data"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  t_ls.push_back(new_lock("BD_lengths_present",true));
  lock_tests["write_LAMMPS_data"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["write_LAMMPS_data"] = t_ls;

  // spherical_bdry
  // number of required parameters
  N_param_reqs["spherical_bdry"] = 4;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("BD_lengths_present",true));
  lock_tests["spherical_bdry"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["spherical_bdry"] = t_ls;

  // switch_bonds
  // number of required parameters
  N_param_reqs["switch_bonds"] = 1;
  // lock tests
  t_ls.clear();
  lock_tests["switch_bonds"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["switch_bonds"] = t_ls;

  // switch_bending_angles
  // number of required parameters
  N_param_reqs["switch_bending_angles"] = 1;
  // lock tests
  t_ls.clear();
  lock_tests["switch_bending_angles"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["switch_bending_angles"] = t_ls;
  
  // switch_twisting_angles
  // number of required parameters
  N_param_reqs["switch_twisting_angles"] = 1;
  // lock tests
  t_ls.clear();
  lock_tests["switch_twisting_angles"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["switch_twisting_angles"] = t_ls;

  // write_mono_xyz
  // number of required parameters
  N_param_reqs["write_mono_xyz"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["write_mono_xyz"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["write_mono_xyz"] = t_ls;

  
  /////////////////////////////////
  // mapper command requirements //
  /////////////////////////////////
  

  // set_initial_state
  // number of required parameters
  N_param_reqs["set_initial_state"] = 0;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["set_initial_state"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("map_initial_present",true));
  lock_updates["set_initial_state"] = t_ls;

  // set_final_state
  // number of required parameters
  N_param_reqs["set_final_state"] = 0;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  t_ls.push_back(new_lock("map_initial_present",true));
  lock_tests["set_final_state"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("map_final_present",true));
  lock_updates["set_final_state"] = t_ls;

  // map_replication
  // number of required parameters
  N_param_reqs["map_replication"] = 0;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("map_initial_present",true));
  t_ls.push_back(new_lock("map_final_present",true));
  t_ls.push_back(new_lock("BD_lengths_present",true));
  lock_tests["map_replication"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["map_replication"] = t_ls;


  ///////////////////////////////////////////
  // LAMMPS_simulator command requirements //
  ///////////////////////////////////////////


  // prepare_simulator
  // number of required parameters
  N_param_reqs["prepare_simulator"] = 1;
  // lock tests
  t_ls.clear();
  lock_tests["prepare_simulator"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_updates["prepare_simulator"] = t_ls;

  // simulator_include_file
  // number of required parameters
  N_param_reqs["simulator_include_file"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["simulator_include_file"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_include_file"] = t_ls;

  // sync_simulator_and_system
  // number of required parameters
  N_param_reqs["sync_simulator_and_system"] = 0;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["sync_simulator_and_system"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["sync_simulator_and_system"] = t_ls;

  // clear_simulator
  // number of required parameters
  N_param_reqs["clear_simulator"] = 0;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["clear_simulator"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["clear_simulator"] = t_ls;

  // simulator_set_nProc
  // number of required parameters
  N_param_reqs["simulator_set_nProc"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["simulator_set_nProc"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_set_nProc"] = t_ls;

  // simulator_set_prng_seed
  // number of required parameters
  N_param_reqs["simulator_set_prng_seed"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["simulator_set_prng_seed"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_set_prng_seed"] = t_ls;

  // simulator_set_DNA_model
  // number of required parameters
  N_param_reqs["simulator_set_DNA_model"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["simulator_set_DNA_model"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("DNA_model",true));
  lock_updates["simulator_set_DNA_model"] = t_ls;

  // simulator_set_output_details
  // number of required parameters
  N_param_reqs["simulator_set_output_details"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["simulator_set_output_details"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("output_details",true));
  lock_updates["simulator_set_output_details"] = t_ls;

  // simulator_set_delta_t
  // number of required parameters
  N_param_reqs["simulator_set_delta_t"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["simulator_set_delta_t"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("delta_t",true));
  lock_updates["simulator_set_delta_t"] = t_ls;

  // simulator_store_timestep
  // number of required parameters
  N_param_reqs["simulator_store_timestep"] = 0;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["simulator_store_timestep"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_store_timestep"] = t_ls;

  // simulator_restore_timestep
  // number of required parameters
  N_param_reqs["simulator_restore_timestep"] = 0;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["simulator_restore_timestep"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_restore_timestep"] = t_ls;

  // simulator_increment_timestep
  // number of required parameters
  N_param_reqs["simulator_increment_timestep"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["simulator_increment_timestep"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_increment_timestep"] = t_ls;

  // simulator_reset_timestep
  // number of required parameters
  N_param_reqs["simulator_reset_timestep"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["simulator_reset_timestep"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_reset_timestep"] = t_ls;

  // simulator_reset_prev_dump_timestep
  // number of required parameters
  N_param_reqs["simulator_reset_prev_dump_timestep"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  lock_tests["simulator_reset_prev_dump_timestep"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_reset_prev_dump_timestep"] = t_ls;

  // simulator_read_data
  // number of required parameters
  N_param_reqs["simulator_read_data"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("DNA_model",true));
  t_ls.push_back(new_lock("output_details",true));
  t_ls.push_back(new_lock("delta_t",true));
  lock_tests["simulator_read_data"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_updates["simulator_read_data"] = t_ls;

  // simulator_minimize_soft_harmonic
  // number of required parameters
  N_param_reqs["simulator_minimize_soft_harmonic"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_tests["simulator_minimize_soft_harmonic"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_minimize_soft_harmonic"] = t_ls;

  // simulator_minimize_hard_harmonic
  // number of required parameters
  N_param_reqs["simulator_minimize_hard_harmonic"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_tests["simulator_minimize_hard_harmonic"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_minimize_hard_harmonic"] = t_ls;

  // simulator_minimize_topoDNA_harmonic
  // number of required parameters
  N_param_reqs["simulator_minimize_topoDNA_harmonic"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_tests["simulator_minimize_topoDNA_harmonic"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_minimize_topoDNA_harmonic"] = t_ls;

  // simulator_minimize_soft_FENE
  // number of required parameters
  N_param_reqs["simulator_minimize_soft_FENE"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_tests["simulator_minimize_soft_FENE"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_minimize_soft_FENE"] = t_ls;

  // simulator_minimize_hard_FENE
  // number of required parameters
  N_param_reqs["simulator_minimize_hard_FENE"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_tests["simulator_minimize_hard_FENE"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_minimize_hard_FENE"] = t_ls;

  // simulator_minimize_topoDNA_FENE
  // number of required parameters
  N_param_reqs["simulator_minimize_topoDNA_FENE"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_tests["simulator_minimize_topoDNA_FENE"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_minimize_topoDNA_FENE"] = t_ls;

  // simulator_run_soft_harmonic
  // number of required parameters
  N_param_reqs["simulator_run_soft_harmonic"] = 5;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_tests["simulator_run_soft_harmonic"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_run_soft_harmonic"] = t_ls;

  // simulator_run_hard_harmonic
  // number of required parameters
  N_param_reqs["simulator_run_hard_harmonic"] = 5;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_tests["simulator_run_hard_harmonic"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_run_hard_harmonic"] = t_ls;

  // simulator_run_topoDNA_harmonic
  // number of required parameters
  N_param_reqs["simulator_run_topoDNA_harmonic"] = 5;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_tests["simulator_run_topoDNA_harmonic"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_run_topoDNA_harmonic"] = t_ls;

  // simulator_run_soft_FENE
  // number of required parameters
  N_param_reqs["simulator_run_soft_FENE"] = 5;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_tests["simulator_run_soft_FENE"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_run_soft_FENE"] = t_ls;

  // simulator_run_hard_FENE
  // number of required parameters
  N_param_reqs["simulator_run_hard_FENE"] = 5;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_tests["simulator_run_hard_FENE"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_run_hard_FENE"] = t_ls;

  // simulator_run_topoDNA_FENE
  // number of required parameters
  N_param_reqs["simulator_run_topoDNA_FENE"] = 5;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_tests["simulator_run_topoDNA_FENE"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_run_topoDNA_FENE"] = t_ls;

  // simulator_load_loop_params
  // number of required parameters
  N_param_reqs["simulator_load_loop_params"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  // t_ls.push_back(new_lock("lmp_data_present",true));
  lock_tests["simulator_load_loop_params"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("loop_params_present",true));
  lock_updates["simulator_load_loop_params"] = t_ls;

  // simulator_run_loops
  // number of required parameters
  N_param_reqs["simulator_run_loops"] = 6;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("lmp_data_present",true));
  t_ls.push_back(new_lock("loop_params_present",true));
  lock_tests["simulator_run_loops"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_run_loops"] = t_ls;

  // switch_Ori_bdry_attraction
  // number of required parameters
  N_param_reqs["switch_Ori_bdry_attraction"] = 1;
  // lock tests
  t_ls.clear();
  lock_tests["switch_Ori_bdry_attraction"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["switch_Ori_bdry_attraction"] = t_ls;

  // switch_Ori_pair_repulsion
  // number of required parameters
  N_param_reqs["switch_Ori_pair_repulsion"] = 1;
  // lock tests
  t_ls.clear();
  lock_tests["switch_Ori_pair_repulsion"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["switch_Ori_pair_repulsion"] = t_ls;
  

  /////////////////////////////////
  // fused command requirements //
  ////////////////////////////////


  // sys_write_sim_read_LAMMPS_data
  // number of required parameters
  N_param_reqs["sys_write_sim_read_LAMMPS_data"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  t_ls.push_back(new_lock("BD_lengths_present",true));
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("DNA_model",true));
  t_ls.push_back(new_lock("output_details",true));
  t_ls.push_back(new_lock("delta_t",true));
  lock_tests["sys_write_sim_read_LAMMPS_data"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_updates["sys_write_sim_read_LAMMPS_data"] = t_ls;

  // sys_write_sim_read_LAMMPS_data_at_timestep
  // number of required parameters
  N_param_reqs["sys_write_sim_read_LAMMPS_data_at_timestep"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  t_ls.push_back(new_lock("BD_lengths_present",true));
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("DNA_model",true));
  t_ls.push_back(new_lock("output_details",true));
  t_ls.push_back(new_lock("delta_t",true));
  lock_tests["sys_write_sim_read_LAMMPS_data_at_timestep"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_updates["sys_write_sim_read_LAMMPS_data_at_timestep"] = t_ls;

  // simulator_relax_progressive
  // number of required parameters
  N_param_reqs["simulator_relax_progressive"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_tests["simulator_relax_progressive"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_relax_progressive"] = t_ls;

  // simulator_expand_bdry_particles
  // number of required parameters
  N_param_reqs["simulator_expand_bdry_particles"] = 4;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("simulator_prepared",true));
  t_ls.push_back(new_lock("lmp_data_present",true));
  lock_tests["simulator_expand_bdry_particles"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["simulator_expand_bdry_particles"] = t_ls;
  

  //////////////////////////////
  // metacommand requirements //
  //////////////////////////////
  
  
  // repeat
  // number of required parameters
  N_param_reqs["repeat"] = 1;
  // lock tests
  t_ls.clear();
  lock_tests["repeat"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["repeat"] = t_ls;

  // end_repeat
  // number of required parameters
  N_param_reqs["end_repeat"] = 0;
  // lock tests
  t_ls.clear();
  lock_tests["end_repeat"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["end_repeat"] = t_ls;

  // repeat_replicates
  // number of required parameters
  N_param_reqs["repeat_replicates"] = 3;
  // lock tests
  t_ls.clear();
  lock_tests["repeat_replicates"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["repeat_replicates"] = t_ls;

  // end_repeat_replicates
  // number of required parameters
  N_param_reqs["end_repeat_replicates"] = 0;
  // lock tests
  t_ls.clear();
  lock_tests["end_repeat_replicates"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["end_repeat_replicates"] = t_ls;
  
}


// function to execute commands
int btree_driver::execute_commands()
{
  
  std::cout << "\n---BEGIN EXECUTING COMMANDS---\n" << std::endl;

  reset_command_locks_and_updates();
  
  int e = 0;
    
  for (size_t i_c=0; i_c<commands.size(); i_c++)
    {
      e += test_command_lock_validity(commands[i_c]);
      if (e > 0) break;
      e += execute_single_command(commands[i_c],command_params[i_c]);
      if (e > 0) break;
      update_lock_state_post_command(commands[i_c]);
    }

  std::cout << "---END EXECUTING COMMANDS---\n" << std::endl;

  return e;
}


// function to execute single command
int btree_driver::execute_single_command(std::string &command,
					 std::vector<std::string> &params)
{

  std::cout << "\nCOMMAND: " << command << std::endl;
  for (size_t i_p=0; i_p<params.size(); i_p++)
    {
      std::cout << "\tparam_" << i_p
	   << ": " << params[i_p]
	   << std::endl;
    }
  std::cout << "\n" << std::endl;


  int error_code = 0;
  
  //////////////////
  // COMMAND LIST //
  //////////////////


  if (command == "terminate")
    {
      error_code = terminate();
    }


  else if (command == "switch_skip_runs")
    {
      error_code = switch_skip_runs(params);
    }      


  //////////////////////////
  // Binary Tree Commands //
  //////////////////////////


  // seed the PRNG for the btree
  else if (command == "btree_prng_seed")
    {
      error_code = btree_prng_seed(params);
    }      

      
  // create a new chromosome
  else if (command == "new_chromo")
    {
      error_code = new_chromo(params);
    }
      
      
  // read input state from a file
  else if (command == "input_state")
    {
      error_code = input_state(params);
    }
      

  // write the state to an output file
  else if (command == "output_state")
    {
      error_code = output_state(params);
    }


  // write the state to an output file
  else if (command == "output_state_at_timestep")
    {
      error_code = output_state_at_timestep(params);
    }

      
  // apply state transformations from a file
  else if (command == "transforms_file")
    {
      error_code = transforms_file(params);
    }

      
  // apply state transformations from a file
  else if (command == "transform")
    {
      error_code = transform(params);
    }

      
  // apply random transformations
  else if (command == "random_transforms")
    {
      error_code = random_transforms(params);
    }

      
  // apply state transformations from a file
  else if (command == "regions_file")
    {
      error_code = regions_file(params);
    }

      
  // write the region counts to an output file
  else if (command == "dump_regions")
    {
      error_code = dump_regions(params);
    }


  // write the region counts to an output file at the current timestep
  else if (command == "dump_regions_at_timestep")
    {
      error_code = dump_regions_at_timestep(params);
    }

      
  // solve the topology of the current state
  else if (command == "update_topology")
    {
      error_code = update_topology();
    }


  // write the topology to an output file
  else if (command == "dump_topology")
    {
      error_code = dump_topology(params);
    }


  // write the topology to an output file at the current timestep
  else if (command == "dump_topology_at_timestep")
    {
      error_code = dump_topology_at_timestep(params);
    }


  // write the topology to an output file
  else if (command == "dump_fork_partitions")
    {
      error_code = dump_fork_partitions(params);
    }


  // write the topology to an output file at the current timestep
  else if (command == "dump_fork_partitions_at_timestep")
    {
      error_code = dump_fork_partitions_at_timestep(params);
    }

      
  // update the CG map
  else if (command == "update_CG_map")
    {
      error_code = update_CG_map(params);
    }

      
  // dump the CG map
  else if (command == "dump_CG_map")
    {
      error_code = dump_CG_map(params);
    }


  // dump the CG map at the current timestep
  else if (command == "dump_CG_map_at_timestep")
    {
      error_code = dump_CG_map_at_timestep(params);
    }


  // print the current state
  else if (command == "print")
    {
      error_code = print_state();
    }
  

  /////////////////////////
  // LAMMPS_sys Commands //
  /////////////////////////
      

  // load file containing the monomer coordinates
  else if (command == "load_mono_coords")
    {
      error_code = load_mono_coords(params);
    }


  // load file containing the monomer quaternions
  else if (command == "load_mono_quats")
    {
      error_code = load_mono_quats(params);
    }


  // load file containing the ribosome coordinates
  else if (command == "load_ribo_coords")
    {
      error_code = load_ribo_coords(params);
    }


  // load file containing the ribosome quaternions
  else if (command == "load_ribo_quats")
    {
      error_code = load_ribo_quats(params);
    }


  // load file containing the boundary coordinates
  else if (command == "load_bdry_coords")
    {
      error_code = load_bdry_coords(params);
    }


  // write file containing the monomer coordinates
  else if (command == "write_mono_coords")
    {
      error_code = write_mono_coords(params);
    }


  // write file containing the monomer quaternions
  else if (command == "write_mono_quats")
    {
      error_code = write_mono_quats(params);
    }


  // write file containing the ribosome coordinates
  else if (command == "write_ribo_coords")
    {
      error_code = write_ribo_coords(params);
    }


  // write file containing the ribosome quaternions
  else if (command == "write_ribo_quats")
    {
      error_code = write_ribo_quats(params);
    }


  // write file containing the boundary coordinates
  else if (command == "write_bdry_coords")
    {
      error_code = write_bdry_coords(params);
    }


  // write file containing the boundary coordinates
  else if (command == "spherical_bdry")
    {
      error_code = spherical_bdry(params);
    }


  // load file containing length-scales for BD simulations
  else if (command == "load_BD_lengths")
    {
      error_code = load_BD_lengths(params);
    }


  // switch bonds on/off
  else if (command == "switch_bonds")
    {
      error_code = switch_bonds(params);
    }


  // switch bending angles on/off
  else if (command == "switch_bending_angles")
    {
      error_code = switch_bending_angles(params);
    }


  // switch twisting angles on/off
  else if (command == "switch_twisting_angles")
    {
      error_code = switch_twisting_angles(params);
    }


  // set initial state for mapper
  else if (command == "set_initial_state")
    {
      error_code = set_initial_state();
    }


  // set final state for mapper
  else if (command == "set_final_state")
    {
      error_code = set_final_state();
    }


  // map the replication based on the initial and final states
  else if (command == "map_replication")
    {
      error_code = map_replication();
    }
      

  // write the LAMMPS system data
  else if (command == "write_LAMMPS_data")
    {
      error_code = write_LAMMPS_data(params);
    }

      
  // write the LAMMPS system data
  else if (command == "write_mono_xyz")
    {
      error_code = write_mono_xyz(params);
    }

      
  ////////////////////////
  // Simulator Commands //
  ////////////////////////
      

  // prepare the simulator
  else if (command == "prepare_simulator")
    {
      error_code = prepare_simulator(params);
    }


  // sync the simulator and the system
  else if (command == "sync_simulator_and_system")
    {
      error_code = sync_simulator_and_system();
    }


  // sync the simulator and the system
  else if (command == "clear_simulator")
    {
      error_code = clear_simulator();
    }


  // run a file using the simulator
  else if (command == "simulator_include_file")
    {
      error_code = simulator_include_file(params);
    }


  // read data (LAMMPS data.* file) into the simulator
  else if (command == "simulator_read_data")
    {
      error_code = simulator_read_data(params);
    }


  // set the prng_seed for the simulator
  else if (command == "simulator_set_prng_seed")
    {
      error_code = simulator_set_prng_seed(params);
    }


  // set the number of processors for the simulator
  else if (command == "simulator_set_nProc")
    {
      error_code = simulator_set_nProc(params);
    }


  // set the DNA model for the simulator
  else if (command == "simulator_set_DNA_model")
    {
      error_code = simulator_set_DNA_model(params);
    }


  // set the output details for the simulator
  else if (command == "simulator_set_output_details")
    {
      error_code = simulator_set_output_details(params);
    }


  // set the timestep for the simulator
  else if (command == "simulator_set_delta_t")
    {
      error_code = simulator_set_delta_t(params);
    }


  // store the simulator's timestep
  else if (command == "simulator_store_timestep")
    {
      error_code = simulator_store_timestep();
    }


  // restore the simulator's timestep from the stored value
  else if (command == "simulator_restore_timestep")
    {
      error_code = simulator_restore_timestep();
    }


  // increment simulator's timestep
  else if (command == "simulator_increment_timestep")
    {
      error_code = simulator_increment_timestep(params);
    }


  // reset the simulator's timestep to the specified value
  else if (command == "simulator_reset_timestep")
    {
      error_code = simulator_reset_timestep(params);
    }


  // reset the simulator's previous dump timestep
  else if (command == "simulator_reset_prev_dump_timestep")
    {
      error_code = simulator_reset_prev_dump_timestep(params);
    }


  // minimize with soft potentials and harmonic bonds
  else if (command == "simulator_minimize_soft_harmonic")
    {
      error_code = simulator_minimize<0,0>(params);
    }


  // minimize with hard potentials and harmonic bonds
  else if (command == "simulator_minimize_hard_harmonic")
    {
      error_code = simulator_minimize<1,0>(params);
    }


  // minimize with topoDNA potentials and harmonic bonds
  else if (command == "simulator_minimize_topoDNA_harmonic")
    {
      error_code = simulator_minimize<2,0>(params);
    }


  // minimize with soft potentials and FENE bonds
  else if (command == "simulator_minimize_soft_FENE")
    {
      error_code = simulator_minimize<0,1>(params);
    }


  // minimize with hard potentials and FENE bonds
  else if (command == "simulator_minimize_hard_FENE")
    {
      error_code = simulator_minimize<1,1>(params);
    }


  // minimize with topoDNA potentials and FENE bonds
  else if (command == "simulator_minimize_topoDNA_FENE")
    {
      error_code = simulator_minimize<2,1>(params);
    }


  // run with soft potentials and harmonic bonds
  else if (command == "simulator_run_soft_harmonic")
    {
      error_code = simulator_run<0,0>(params);
    }


  // run with hard potentials and harmonic bonds
  else if (command == "simulator_run_hard_harmonic")
    {
      error_code = simulator_run<1,0>(params);
    }


  // run with topoDNA potentials and harmonic bonds
  else if (command == "simulator_run_topoDNA_harmonic")
    {
      error_code = simulator_run<2,0>(params);
    }


  // run with soft potentials and FENE bonds
  else if (command == "simulator_run_soft_FENE")
    {
      error_code = simulator_run<0,1>(params);
    }


  // run with hard potentials and FENE bonds
  else if (command == "simulator_run_hard_FENE")
    {
      error_code = simulator_run<1,1>(params);
    }


  // run with topoDNA potentials and FENE bonds
  else if (command == "simulator_run_topoDNA_FENE")
    {
      error_code = simulator_run<2,1>(params);
    }

      
  // load file containing loop parameters
  else if (command == "simulator_load_loop_params")
    {
      error_code = simulator_load_loop_params(params);
    }


  // run a simulation with loops
  else if (command == "simulator_run_loops")
    {
      error_code = simulator_run_loops(params);
    }


  // run a simulation with bdry particles that increase in size
  else if (command == "simulator_expand_bdry_particles")
    {
      error_code = simulator_expand_bdry_particles(params);
    }


  // switch on the Ori bdry attraction
  else if (command == "switch_Ori_bdry_attraction")
    {
      error_code = switch_Ori_bdry_attraction(params);
    }


  // switch on the Ori pair_repulsion
  else if (command == "switch_Ori_pair_repulsion")
    {
      error_code = switch_Ori_pair_repulsion(params);
    }


  ////////////////////
  // Fused Commands //
  ////////////////////


  // write LAMMPS data with system and read LAMMPS data with simulator
  else if (command == "sys_write_sim_read_LAMMPS_data")
    {
      error_code = sys_write_sim_read_LAMMPS_data(params);
    }


  // write LAMMPS data with system and read LAMMPS data with simulator
  else if (command == "sys_write_sim_read_LAMMPS_data_at_timestep")
    {
      error_code = sys_write_sim_read_LAMMPS_data_at_timestep(params);
    }

  
  // write LAMMPS data with system and read LAMMPS data with simulator
  else if (command == "simulator_relax_progressive")
    {
      error_code = simulator_relax_progressive(params);
    }


  // the command is unrecognized
  else
    {
      std::cout << "command is unrecognized!" << std::endl;
      error_code = 1;
    }

  return error_code;
  
}


//////////////////////////////////////////
// set of functions to perform commands //
//////////////////////////////////////////


int btree_driver::terminate()
{
  return 1;
}


int btree_driver::switch_skip_runs(std::vector<std::string> &params)
{
  if (params[0] == "T")
    {
      skip_runs = true;
    }
  else if (params[0] == "F")
    {
      skip_runs = false;
    }
  else
    {
      std::cout << "ERROR: invalid switch" << std::endl;
      return 1;
    }
  return 0;
}


int btree_driver::new_chromo(std::vector<std::string> &params)
{
  btree_state driver_st;
  driver_st.size = stoi(params[0]);
  driver_st.transforms.clear();
  driver_bt.prepare_state(driver_st);
  return 0;
}


int btree_driver::input_state(std::vector<std::string> &params)
{
  btree_state driver_st;
  driver_st = driver_bt.read_state(params[0]);
  driver_bt.prepare_state(driver_st);
  return 0;
}


int btree_driver::output_state(std::vector<std::string> &params)
{
  driver_bt.write_state(params[0],driver_bt.get_state());
  return 0;
}


int btree_driver::output_state_at_timestep(std::vector<std::string> &params)
{
  int e;
  std::string ts_mod = get_timestep_modifier();
  std::vector<std::string> params_w_ts;

  for (std::string param : params)
    {
      params_w_ts.push_back(param);
    }

  insert_timestep_modifier(ts_mod,params_w_ts[0]);
  
  e = output_state(params_w_ts);
  
  return e;
}


int btree_driver::transforms_file(std::vector<std::string> &params)
{
  btree_transforms driver_tr;
  driver_tr = driver_bt.read_transforms(params[0]);
  driver_bt.apply_transforms(driver_tr);
  return 0;
}


int btree_driver::transform(std::vector<std::string> &params)
{
  driver_bt.single_transform(driver_bt.parse_transform(params[0]));
  return 0;
}


int btree_driver::random_transforms(std::vector<std::string> &params)
{
  driver_bt.random_transforms(stoi(params[0]));
  return 0;
}


int btree_driver::regions_file(std::vector<std::string> &params)
{
  driver_rg.clear();
  driver_rg = driver_bt.read_regions(params[0],stoi(params[1]));
  return 0;
}


int btree_driver::dump_regions(std::vector<std::string> &params)
{
  // update topology before updating regions
  if (lock_state["topo_update"] == true)
    {
      driver_bt.solve_topology();
      lock_state["topo_update"] = false;
    }
  driver_bt.update_region_counts(driver_rg);
  driver_bt.dump_regions(params[0],driver_rg);
  return 0;
}


int btree_driver::dump_regions_at_timestep(std::vector<std::string> &params)
{
  int e;
  std::string ts_mod = get_timestep_modifier();
  std::vector<std::string> params_w_ts;

  for (std::string param : params)
    {
      params_w_ts.push_back(param);
    }

  insert_timestep_modifier(ts_mod,params_w_ts[0]);
  
  e = dump_regions(params_w_ts);
  
  return e;
}


int btree_driver::dump_topology(std::vector<std::string> &params)
{
  // update topology before dumping
  if (lock_state["topo_update"] == true)
    {
      driver_bt.solve_topology();
      lock_state["topo_update"] = false;
    }
  driver_bt.dump_topology(params[0],stoi(params[1]));
  return 0;
}


int btree_driver::dump_topology_at_timestep(std::vector<std::string> &params)
{
  int e;
  std::string ts_mod = get_timestep_modifier();
  std::vector<std::string> params_w_ts;

  for (std::string param : params)
    {
      params_w_ts.push_back(param);
    }

  insert_timestep_modifier(ts_mod,params_w_ts[0]);
  
  e = dump_topology(params_w_ts);
  
  return e;
}


int btree_driver::dump_fork_partitions(std::vector<std::string> &params)
{
  // update topology before dumping
  if (lock_state["topo_update"] == true)
    {
      driver_bt.solve_topology();
      lock_state["topo_update"] = false;
    }
  driver_bt.dump_fork_partitions(params[0],stoi(params[1]));
  return 0;
}


int btree_driver::dump_fork_partitions_at_timestep(std::vector<std::string> &params)
{
  int e;
  std::string ts_mod = get_timestep_modifier();
  std::vector<std::string> params_w_ts;

  for (std::string param : params)
    {
      params_w_ts.push_back(param);
    }

  insert_timestep_modifier(ts_mod,params_w_ts[0]);
  
  e = dump_fork_partitions(params_w_ts);
  
  return e;
}


int btree_driver::update_topology()
{
  driver_bt.solve_topology();
  // disable flag for topology updating after an update
  lock_state["topo_update"] = false;
  return 0;
}


int btree_driver::update_CG_map(std::vector<std::string> &params)
{
  // update topology before dumping
  if (lock_state["topo_update"] == true)
    {
      driver_bt.solve_topology();
      lock_state["topo_update"] = false;
    }
  driver_CG = driver_bt.update_CG_map(stoi(params[0]));
  return 0;
}


int btree_driver::dump_CG_map(std::vector<std::string> &params)
{
  // update topology before updating CG_map
  if (lock_state["topo_update"] == true)
    {
      driver_bt.solve_topology();
      lock_state["topo_update"] = false;
    }
  // update CG_map before dumping
  if (lock_state["CG_update"] == true)
    {
      driver_CG = driver_bt.update_CG_map(stoi(params[1]));
      lock_state["CG_update"] = false;
    }
  driver_bt.dump_CG_map(params[0],stoi(params[2]),driver_CG);
  return 0;
}


int btree_driver::dump_CG_map_at_timestep(std::vector<std::string> &params)
{
  int e;
  std::string ts_mod = get_timestep_modifier();
  std::vector<std::string> params_w_ts;

  for (std::string param : params)
    {
      params_w_ts.push_back(param);
    }

  insert_timestep_modifier(ts_mod,params_w_ts[0]);
  
  e = dump_CG_map(params_w_ts);
  
  return e;
}


int btree_driver::load_BD_lengths(std::vector<std::string> &params)
{
  driver_lmp_sys.read_BD_lengths(params[0]);
  return 0;
}


int btree_driver::load_mono_coords(std::vector<std::string> &params)
{
  driver_lmp_sys.set_btree(driver_bt.get_state());
  int e = driver_lmp_sys.read_mono_coords(params[0],params[1]);
  return e;
}


int btree_driver::load_mono_quats(std::vector<std::string> &params)
{
  driver_lmp_sys.set_btree(driver_bt.get_state());
  int e = driver_lmp_sys.read_mono_quats(params[0],params[1]);
  return e;
}


int btree_driver::load_ribo_coords(std::vector<std::string> &params)
{
  int e = driver_lmp_sys.read_ribo_coords(params[0],params[1]);
  return e;
}


int btree_driver::load_ribo_quats(std::vector<std::string> &params)
{
  int e = driver_lmp_sys.read_ribo_quats(params[0],params[1]);
  return e;
}


int btree_driver::load_bdry_coords(std::vector<std::string> &params)
{
  int e = driver_lmp_sys.read_bdry_coords(params[0],params[1]);
  return e;
}


int btree_driver::write_mono_coords(std::vector<std::string> &params)
{
  int e = driver_lmp_sys.write_mono_coords(params[0],params[1]);
  return e;
}


int btree_driver::write_mono_quats(std::vector<std::string> &params)
{
  int e = driver_lmp_sys.write_mono_quats(params[0],params[1]);
  return e;
}


int btree_driver::write_ribo_coords(std::vector<std::string> &params)
{
  int e = driver_lmp_sys.write_ribo_coords(params[0],params[1]);
  return e;
}


int btree_driver::write_ribo_quats(std::vector<std::string> &params)
{
  int e = driver_lmp_sys.write_ribo_quats(params[0],params[1]);
  return e;
}


int btree_driver::write_bdry_coords(std::vector<std::string> &params)
{
  int e = driver_lmp_sys.write_bdry_coords(params[0],params[1]);
  return e;
}


int btree_driver::write_LAMMPS_data(std::vector<std::string> &params)
{
  driver_lmp_sys.set_btree(driver_bt.get_state());
  driver_lmp_sys.write_data(params[0]);
  return 0;
}


int btree_driver::spherical_bdry(std::vector<std::string> &params)
{
  driver_lmp_sys.generate_spherical_bdry(stod(params[0]),
					 stod(params[1]),
					 stod(params[2]),
					 stod(params[3]));
  return 0;
}


int btree_driver::switch_bonds(std::vector<std::string> &params)
{
  if (params[0] == "T")
    {
      driver_lmp_sys.switch_bonds(true);
    }
  else if (params[0] == "F")
    {
      driver_lmp_sys.switch_bonds(false);
    }
  else
    {
      std::cout << "ERROR: invalid switch" << std::endl;
      return 1;
    }
  return 0;
}


int btree_driver::switch_bending_angles(std::vector<std::string> &params)
{
  if (params[0] == "T")
    {
      driver_lmp_sys.switch_bending_angles(true);
    }
  else if (params[0] == "F")
    {
      driver_lmp_sys.switch_bending_angles(false);
    }
  else
    {
      std::cout << "ERROR: invalid switch" << std::endl;
      return 1;
    }
  return 0;
}


int btree_driver::switch_twisting_angles(std::vector<std::string> &params)
{
  if (params[0] == "T")
    {
      driver_lmp_sys.switch_twisting_angles(true);
    }
  else if (params[0] == "F")
    {
      driver_lmp_sys.switch_twisting_angles(false);
    }
  else
    {
      std::cout << "ERROR: invalid switch" << std::endl;
      return 1;
    }
  return 0;
}


int btree_driver::switch_extra_potential(std::string extra_pot, std::string s)
{
  if (s == "T")
    {
      driver_lmp_simulator.switch_extra_potential(extra_pot,true);
    }
  else if (s == "F")
    {
      driver_lmp_simulator.switch_extra_potential(extra_pot,false);
    }
  else
    {
      std::cout << "ERROR: invalid switch" << std::endl;
      return 1;
    }
  return 0;
}


int btree_driver::switch_Ori_bdry_attraction(std::vector<std::string> &params)
{
  int e = switch_extra_potential("Ori_bdry_attraction",params[0]);
  return e;
}


int btree_driver::switch_Ori_pair_repulsion(std::vector<std::string> &params)
{
  int e = switch_extra_potential("Ori_pair_repulsion",params[0]);
  return e;
}


int btree_driver::write_mono_xyz(std::vector<std::string> &params)
{
  driver_lmp_sys.write_mono_xyz(params[0]);
  return 0;
}


int btree_driver::set_initial_state()
{
  driver_mapper.set_initial_state(driver_bt.get_state());
  return 0;
}


int btree_driver::set_final_state()
{
  driver_mapper.set_final_state(driver_bt.get_state());
  return 0;
}


int btree_driver::map_replication()
{
  int e = driver_mapper.prepare_mapping();
  driver_lmp_sys.apply_mono_mapping(driver_mapper.get_map());
  return e;
}


int btree_driver::prepare_simulator(std::vector<std::string> &params)
{
  driver_lmp_simulator.LAMMPS_initialize(params[0]);
  driver_lmp_simulator.set_lmp_sys(&driver_lmp_sys);
  return 0;
}


int btree_driver::simulator_include_file(std::vector<std::string> &params)
{
  driver_lmp_simulator.include_file(params[0]);
  return 0;
}


int btree_driver::sync_simulator_and_system()
{
  driver_lmp_sys.set_btree(driver_bt.get_state());
  driver_lmp_simulator.sim_to_sys();
  return 0;
}


int btree_driver::clear_simulator()
{
  driver_lmp_simulator.clear();
  return 0;
}


int btree_driver::simulator_set_nProc(std::vector<std::string> &params)
{
  driver_lmp_simulator.set_nProc(stoi(params[0]));  
  return 0;
}


int btree_driver::simulator_set_prng_seed(std::vector<std::string> &params)
{
  driver_lmp_simulator.set_prng_seed(stoi(params[0]));
  driver_lmp_sys.prng_seed(stoi(params[0]));
  return 0;
}


int btree_driver::simulator_set_DNA_model(std::vector<std::string> &params)
{
  driver_lmp_simulator.set_DNA_model_dir(params[0]);
  return 0;
}


int btree_driver::simulator_set_output_details(std::vector<std::string> &params)
{
  driver_lmp_simulator.set_output_details(params[0],params[1]);
  return 0;
}


int btree_driver::simulator_set_delta_t(std::vector<std::string> &params)
{
  driver_lmp_simulator.set_delta_t(stod(params[0]));
  return 0;
}


int btree_driver::simulator_read_data(std::vector<std::string> &params)
{
  driver_lmp_simulator.clear();
  driver_lmp_simulator.reset_protocol_variables();
  driver_lmp_simulator.global_setup();
  driver_lmp_simulator.read_data(params[0]);
  driver_lmp_simulator.standard_computes();
  return 0;  
}


template<int SOFT_HARD_TOPO, int HARMONIC_FENE>
int btree_driver::simulator_minimize(std::vector<std::string> &params)
{
  thermo_dump_parameters t_d_p;

  t_d_p.append = false;
  t_d_p.write_first = true;
  t_d_p.dump_freq = 0;
  t_d_p.thermo_freq = stoi(params[0]);

  // skip the minimization if skipping runs
  if (skip_runs == true) return 0;

  // run the minimization
  if (HARMONIC_FENE == 0)
    {
      if (SOFT_HARD_TOPO == 0)
	{
	  driver_lmp_simulator.minimize_soft_harmonic(t_d_p);
	}
      else if (SOFT_HARD_TOPO == 1)
	{
	  driver_lmp_simulator.minimize_hard_harmonic(t_d_p);
	}
      else if (SOFT_HARD_TOPO == 2)
	{
	  driver_lmp_simulator.minimize_topoDNA_harmonic(t_d_p);
	}
    }
  else if (HARMONIC_FENE == 1)
    {
      if (SOFT_HARD_TOPO == 0)
	{
	  driver_lmp_simulator.minimize_soft_FENE(t_d_p);
	}
      else if (SOFT_HARD_TOPO == 1)
	{
	  driver_lmp_simulator.minimize_hard_FENE(t_d_p);
	}
      else if (SOFT_HARD_TOPO == 2)
	{
	  driver_lmp_simulator.minimize_topoDNA_FENE(t_d_p);
	}
    }


  return 0;
}


template<int SOFT_HARD_TOPO, int HARMONIC_FENE>
int btree_driver::simulator_run(std::vector<std::string> &params)
{
  thermo_dump_parameters t_d_p;

  t_d_p.append = false;
  t_d_p.write_first = true;
  if (params[3] == "append")
    {
      t_d_p.append = true;
      t_d_p.write_first = false;
    }
  if (params[4] == "skip_first") t_d_p.write_first = false;
  t_d_p.dump_freq = stoi(params[2]);
  t_d_p.thermo_freq = stoi(params[1]);

  // skip the run if skipping runs
  if (skip_runs == true)
    {
      // still increment timestep for testing
      driver_lmp_simulator.increment_Nt(stoul(params[0]));
      return 0;
    }

  // run the Brownian dynamics
  if (HARMONIC_FENE == 0)
    {
      if (SOFT_HARD_TOPO == 0)
	{
	  driver_lmp_simulator.run_soft_harmonic(stoul(params[0]),t_d_p);
	}
      else if (SOFT_HARD_TOPO == 1)
	{
	  driver_lmp_simulator.run_hard_harmonic(stoul(params[0]),t_d_p);
	}
      else if (SOFT_HARD_TOPO == 2)
	{
	  driver_lmp_simulator.run_topoDNA_harmonic(stoul(params[0]),t_d_p);
	}
    }
  else if (HARMONIC_FENE == 1)
    {
      if (SOFT_HARD_TOPO == 0)
	{
	  driver_lmp_simulator.run_soft_FENE(stoul(params[0]),t_d_p);
	}
      else if (SOFT_HARD_TOPO == 1)
	{
	  driver_lmp_simulator.run_hard_FENE(stoul(params[0]),t_d_p);
	}
      else if (SOFT_HARD_TOPO == 2)
	{
	  driver_lmp_simulator.run_topoDNA_FENE(stoul(params[0]),t_d_p);
	}
    }

  return 0;
}


int btree_driver::simulator_load_loop_params(std::vector<std::string> &params)
{
  int e = driver_lmp_simulator.read_loop_params(params[0]);
  return e;
}


int btree_driver::simulator_run_loops(std::vector<std::string> &params)
{
  thermo_dump_parameters t_d_p;

  t_d_p.append = false;
  t_d_p.write_first = true;
  if (params[4] == "append")
    {
      t_d_p.append = true;
      t_d_p.write_first = false;
    }
  if (params[5] == "skip_first") t_d_p.write_first = false;
  t_d_p.dump_freq = stoi(params[3]);
  t_d_p.thermo_freq = stoi(params[2]);

  // skip the run if skipping runs
  if (skip_runs == true)
    {
      // still increment timestep for testing
      driver_lmp_simulator.increment_Nt(stoul(params[1]));
      return 0;
    }

  driver_lmp_simulator.run_loops(stoi(params[0]),
				 stoul(params[1]),
				 t_d_p);

  return 0;
}


int btree_driver::simulator_store_timestep()
{
  driver_lmp_simulator.store_Nt();
  return 0;
}


int btree_driver::simulator_restore_timestep()
{
  driver_lmp_simulator.restore_Nt();
  return 0;
}


int btree_driver::simulator_increment_timestep(std::vector<std::string> &params)
{
  driver_lmp_simulator.increment_Nt(stoul(params[0]));
  return 0;
}


int btree_driver::simulator_reset_timestep(std::vector<std::string> &params)
{
  driver_lmp_simulator.reset_Nt(stoul(params[0]));
  return 0;
}


int btree_driver::simulator_reset_prev_dump_timestep(std::vector<std::string> &params)
{
  driver_lmp_simulator.reset_prev_dump_Nt(stoul(params[0]));
  return 0;
}


int btree_driver::btree_prng_seed(std::vector<std::string> &params)
{
  // seed the PRNG
  driver_bt.prng_seed(stoi(params[0]));
  return 0;
}


int btree_driver::print_state()
{
  // update topology before printing
  if (lock_state["topo_update"] == true)
    {
      driver_bt.solve_topology();
      lock_state["topo_update"] = false;
    }
  driver_bt.print_tree();
  return 0;
}


int btree_driver::sys_write_sim_read_LAMMPS_data(std::vector<std::string> &params)
{
  int e = 0;
  e += write_LAMMPS_data(params);
  e += simulator_read_data(params);
  return e;
}


int btree_driver::sys_write_sim_read_LAMMPS_data_at_timestep(std::vector<std::string> &params)
{
  int e = 0;
  std::string ts_mod = get_timestep_modifier();
  std::vector<std::string> params_w_ts;

  for (std::string param : params)
    {
      params_w_ts.push_back(param);
    }

  append_timestep_modifier(ts_mod,params_w_ts[0]);
  
  e += write_LAMMPS_data(params_w_ts);
  e += simulator_read_data(params_w_ts);
  return e;
}


int btree_driver::simulator_relax_progressive(std::vector<std::string> &params)
{
  int e = 0;

  std::vector<std::string> min_params;
  std::vector<std::string> run_params;

  int run_steps = stoi(params[0]);
  int thermo_freq = stoi(params[1]);

  // create the parameter vector for minimizations
  min_params.push_back(std::to_string(thermo_freq));
  // create the parameter vector for runs
  run_params.push_back(std::to_string(run_steps));
  run_params.push_back(std::to_string(thermo_freq));
  run_params.push_back("0");
  run_params.push_back("noappend");
  run_params.push_back("first");

  // skip the run if skipping runs
  if (skip_runs == true) return 0;
  
  // store the current timestep
  e += simulator_store_timestep();
  // minimize with soft pairs and harmonic bonds
  e += simulator_minimize<0,0>(min_params);
  // run with soft pairs and harmonic bonds
  if (run_steps > 0)
    {
      e += simulator_run<0,0>(run_params);
    }
  // minimize with hard pairs and harmonic bonds
  e += simulator_minimize<1,0>(min_params);
  // run with hard pairs and harmonic bonds
  if (run_steps > 0)
    {
      e += simulator_run<1,0>(run_params);
    }
  // minimize with hard pairs and FENE bonds
  e += simulator_minimize<1,1>(min_params);
  // restore the timestep
  e += simulator_restore_timestep();

  return e;
}


int btree_driver::simulator_expand_bdry_particles(std::vector<std::string> &params)
{
  int e = 0;
  double ds;

  std::vector<std::string> run_params;

  double expansion_scale = stod(params[0]);
  int N_iter = stoi(params[1]);
  int run_steps = stoi(params[2]);
  int thermo_freq = stoi(params[3]);

  // create the parameter vector for runs
  run_params.push_back(std::to_string(run_steps));
  run_params.push_back(std::to_string(thermo_freq));
  run_params.push_back("0");
  run_params.push_back("noappend");
  run_params.push_back("first");


  std::cout << "expanding bdry particles by factor of "
	    << expansion_scale << std::endl;
  std::cout << "N_iter = "
	    << N_iter << std::endl;

  // store the current timestep
  e += simulator_store_timestep();

  // iterate over the expansion steps
  for (int i_step=0; i_step<N_iter; i_step++)
    {

      // set the expansion step size
      ds = ((expansion_scale-1.0)/N_iter)*(i_step+1);

      // expand the bdry particles
      driver_lmp_simulator.expand_bdry_particles(ds);
      
      // run the system
      if (run_steps > 0)
	{
	  e += simulator_run<1,1>(run_params);
	}

    }

  // reset bdry particle expansion
  driver_lmp_simulator.reset_bdry_particle_expansion();

  // restore the timestep
  e += simulator_restore_timestep();

  return e;
}
