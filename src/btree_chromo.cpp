#include <btree_chromo.hpp>

int main(int argc, char *argv[])
{

  int error_code = 0;
  btree_driver driver;

  if (argc != 2)
    {
      std::cout << "ERROR: missing directive file as command-line argument" << std::endl;
      return 1;
    }

  std::string drctv_filename(argv[argc-1]);

  // read the directives from the directive file
  driver.read_directives(drctv_filename);

  // print the directives to be executed
  driver.print_directives();

  // parse the directives into commands and parameters
  driver.parse_directives();

  // validate the numbers of parameters
  error_code = driver.validate_command_sequence_parameters();
  if (error_code != 0) return 1;
  
  // expand the metacommands
  error_code = driver.expand_metacommands();
  if (error_code != 0) return 1;

  // validate the command sequence
  error_code = driver.validate_command_sequence();
  if (error_code != 0) return 1;

  // print the set of commands
  driver.print_commands();
    
  // execute the commands
  error_code = driver.execute_commands();
  if (error_code != 0)
    {
      std::cout << "error during command execution" << std::endl;
      return 1;
    } 

  return 0;
  
}
