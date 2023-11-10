#include <btree.hpp>

// constructor
btree::btree()
{
  rand_eng.seed(0);
  reset_root();
}

void btree::reset_root()
{
  root = nullptr;
}

// destructor
btree::~btree()
{
  // std::cout << "btree destructor before destroy" << std::endl;
  destroy_tree();
  // std::cout << "btree destructor after destroy" << std::endl;
}


void btree::prng_seed(int s)
{
  rand_eng.seed(s);
}


// function to initialize binary tree
void btree::initialize_tree(int s)
{
  if (root != nullptr) destroy_tree();
  root = new node;
  initialize_branch(root,0,s);
}


// function to prepare a state of the binary tree
void btree::prepare_state(btree_state st)
{

  int error_code;
  
  initialize_tree(st.size);
  for (fork_rho f_r: st.transforms)
    {
      
      error_code = grow_at_branch_asym(f_r.fork,f_r.rho_cw,f_r.rho_ccw);

      if (error_code == -1)
	{
	  std::cout << "ERROR: impossible state, destroying tree" << std::endl;
	  destroy_tree();
	  return;
	}
      
    }
}


// function to get the state of the binary tree
btree_state btree::get_state()
{

  btree_state st;
  fork_rho f_r;
  node *temp_branch;
  
  if (root != nullptr)
    {
      st.size = root->size;

      // add completed forks to state
      for (std::string s: get_completed_forks())
	{
	  
	  f_r.fork = s;
	  temp_branch = get_branch(s);
	  f_r.rho_cw = temp_branch->rho_cw;
	  f_r.rho_ccw = temp_branch->rho_ccw;

	  st.transforms.push_back(f_r);
	  
	}

      // add active forks to state
      for (std::string s: get_active_forks())
	{
	  
	  f_r.fork = s;
	  temp_branch = get_branch(s);
	  f_r.rho_cw = temp_branch->rho_cw;
	  f_r.rho_ccw = temp_branch->rho_ccw;

	  st.transforms.push_back(f_r);
	  
	}
      
    }
  else
    {
      st.size = -1;
    }
  
  return st;
}


// function to apply transforms to the binary tree
void btree::apply_transforms(btree_transforms tr)
{

  int error_code;
  btree_state st = get_state();
  
  for (fork_rho f_r: tr)
    {
      
      error_code = grow_at_branch_asym(f_r.fork,f_r.rho_cw,f_r.rho_ccw);

      if (error_code == -1)
	{
	  std::cout << "ERROR: impossible state, destroying tree and returning to initial state" << std::endl;

	  // destroy the impossible tree and return to the initial state
	  destroy_tree();
	  prepare_state(st);
	  return;
	}
      
    }
}



// function to apply a single transform to the binary tree
void btree::single_transform(fork_rho f_r)
{

  int error_code;
  btree_state st = get_state();
  
  error_code = grow_at_branch_asym(f_r.fork,f_r.rho_cw,f_r.rho_ccw);

  if (error_code == -1)
    {
      std::cout << "ERROR: impossible state, destroying tree and returning to initial state" << std::endl;

      // destroy the impossible tree and return to the initial state
      destroy_tree();
      prepare_state(st);
      return;
    }
      
}


// function to write a binary tree state to a file
void btree::write_state(std::string st_filename, btree_state st)
{
  std::fstream st_file;

  st_file.open(st_filename, std::ios::out);

  if (!st_file.is_open())
    {
      std::cout << "ERROR: file not opened in write_state" << std::endl;
    }
  else
    {

      st_file << "size=" << st.size << std::endl;

      for (fork_rho f_r: st.transforms)
	{
      
	  st_file << f_r.fork <<
	    "_cw" << f_r.rho_cw <<
	    "_ccw" << f_r.rho_ccw <<
	    std::endl;
      
	}
      st_file.close();
    }
}


// function to read a binary tree state from a file
btree_state btree::read_state(std::string st_filename)
{
  std::fstream st_file;

  btree_state st;
  fork_rho f_r;
  
  std::string line;

  int delim;
  
  bool size_found = false;

  st_file.open(st_filename, std::ios::in);

  // std::cout << st_filename << std::endl;

  if (!st_file.is_open())
    {
      st.size = -1;
      std::cout << "ERROR: file not opened in read_state" << std::endl;
    }
  else
    {
      while (1)
	{
	  st_file >> line;
	  if (st_file.eof()) break;
	  
	  // std::cout << line << std::endl;

	  if (size_found == false)
	    {
	      
	      delim = line.find("=");
	      // std::cout << line.substr(0,delim) << std::endl;
	      // std::cout << line.substr(delim+1,line.length()) << std::endl;
	      
	      if (line.substr(0,delim) == "size")
		{
		  st.size = stoi(line.substr(delim+1,line.length()));
		  size_found = true;
		}
	      
	    }
	  else
	    {
	      
	      f_r = parse_transform(line);
	      st.transforms.push_back(f_r);
	      
	    }
      
	}

      st_file.close();
    }

  return st;
}


// function to read a binary tree state from a file
btree_transforms btree::read_transforms(std::string tr_filename)
{
  std::fstream tr_file;

  btree_transforms tr;
  fork_rho f_r;
  
  std::string line;

  tr_file.open(tr_filename, std::ios::in);

  // std::cout << tr_filename << std::endl;

  if (!tr_file.is_open())
    {
      std::cout << "ERROR: file not opened in read_transforms" << std::endl;
    }
  else
    {
      while (1)
	{
	  tr_file >> line;
	  if (tr_file.eof()) break;
	  
	  // std::cout << line << std::endl;
	  f_r = parse_transform(line);
	  tr.push_back(f_r);

     	}

      tr_file.close();
      
    }

  return tr;
}



// function to parse std::string into the transform
fork_rho btree::parse_transform(std::string s)
{
  fork_rho f_r;
  int delim;

  delim = s.find("_");
  if (delim != -1)
    {
      f_r.fork = s.substr(0,delim);
      s.erase(0,delim+1);
    }
  delim = s.find("_");
  if (delim != -1)
    {
      f_r.rho_cw = stoi(s.substr(2,delim));
      s.erase(0,delim+1);
    }
  
  f_r.rho_ccw = stoi(s.substr(3,s.length()));

  return f_r;
  
}


// function to initialize a branch
void btree::initialize_branch(node *branch, int g, int s)
{
  branch->gen = g;
  branch->size = s;
  branch->rho_t = 0;
  branch->rho_cw = 0;
  branch->rho_ccw = 0;
  branch->leaf = true;
  branch->complete = false;
  branch->parent = nullptr;
  branch->left = nullptr;
  branch->right = nullptr;
  branch->topo.start = -1;
  branch->topo.end = -1;
  branch->topo.mid = -1;
  branch->topo.start_link = -1;
  branch->topo.end_link = -1;
}


// function to begin a branching procedure if it is valid
int btree::branch(std::string loc)
{
  node *branch;

  branch = get_branch(loc);

  if (branch == nullptr)
    {
      std::cout << "invalid branch ("
	   << loc
	   << ")"
	   << std::endl;
      
      return 1;
    }
  else
    {
      
      // std::cout << "branching at ("
      // 	   << loc
      // 	   << ")"
      // 	   << std::endl;
      
      split_branch(branch);
      
      return 0;
    }
}


// function to create children at branch
void btree::split_branch(node *branch)
{
  if (branch->leaf == true)
    {
      branch->leaf = false;
  
      branch->left = new node;
      initialize_branch(branch->left,branch->gen+1,branch->size);
      branch->left->parent = branch;
  
      branch->right = new node;
      initialize_branch(branch->right,branch->gen+1,branch->size);
      branch->right->parent = branch;
    }
}


// function to increase rho in a manner dependent on the possible growth
int btree::grow_at_branch_sym(std::string loc, int r_sym)
{
  int rem, error_code;
  int r_cw, r_ccw;

  error_code = branch(loc);
  if (error_code < 0)
    {
      return error_code;
    }

  // grow the branch if r > 0
  if (r_sym > 0)
    {

      r_cw = r_sym/2;
      r_ccw = r_sym - r_cw;

      std::uniform_real_distribution<double> u_dist(0.0,1.0);
      double ur = u_dist(rand_eng);
      if (ur < 0.5)
	{
	  r_ccw = r_cw;
	  r_cw = r_sym - r_ccw;
	}

      rem = grow_at_branch_asym(loc,r_cw,r_ccw);
      
    }
  else // do nothing if r == 0
    {
      rem = 0;
    }
  
  return rem; // return the remainder of the growth
}



// function to increase rho in a manner dependent on the possible growth
int btree::grow_at_branch_asym(std::string loc, int r_cw, int r_ccw)
{
  int error_code;
  node *g_branch;
  std::array<int,2> part_growths {0, 0};

  error_code = branch(loc);
  if (error_code < 0)
    {
      return error_code;
    }

  // grow the branch if r_cw > 0 or r_ccw > 0
  if ((r_cw > 0) || (r_ccw > 0))
    {
      g_branch = get_branch(loc);

      part_growths = partition_growths_sym(g_branch,r_cw,r_ccw);

      g_branch->rho_cw += part_growths[0];
      g_branch->rho_ccw += part_growths[1];
      g_branch->rho_t = g_branch->rho_cw + g_branch->rho_ccw;

      if (g_branch->rho_t >= g_branch->size)
	{
	  g_branch->complete = true;
	  g_branch->rho_cw = g_branch->size;
	  g_branch->rho_ccw = g_branch->size;
	  g_branch->rho_t = g_branch->size;
	}
    
    }
  
  return r_cw + r_ccw - part_growths[0] - part_growths[1]; // return the remainder of the growth
}


// function to calculate maximum growth in clockwise direction
int btree::get_max_growth_cw(node *branch)
{
  int max_size_cw;
  
  if (branch->parent != nullptr)
    {
      if (branch->parent->complete == false)
	{
	  max_size_cw = branch->parent->rho_cw - 1;
	}
      else
	{
	  max_size_cw = branch->size;
	}
    }
  else
    {
      max_size_cw = branch->size - branch->rho_ccw;
    }

  return max_size_cw - branch->rho_cw;
}

// function to calculate maximum growth in counter-clockwise direction
int btree::get_max_growth_ccw(node *branch)
{
  int max_size_ccw;
  
  if (branch->parent != nullptr)
    {
      if (branch->parent->complete == false)
	{
	  max_size_ccw = branch->parent->rho_ccw - 1;
	}
      else
	{
	  max_size_ccw = branch->size;
	}
    }
  else
    {
      max_size_ccw = branch->size - branch->rho_cw;
    }

  return max_size_ccw - branch->rho_ccw;
}

// function to partition maximum growths assuming symmetric rates along cw and ccw
std::array<int,2> btree::partition_growths_sym(node *branch, int proposed_r_cw, int proposed_r_ccw)
{
  std::array<int,2> growths;
  int max_growth_cw, max_growth_ccw;
  
  bool cw_first = true;

  growths[0] = 0;
  growths[1] = 0;

  max_growth_cw = get_max_growth_cw(branch);
  max_growth_ccw = get_max_growth_ccw(branch);

  if (branch->rho_cw < branch->rho_ccw)
    {
      cw_first = true;
    }
  else if (branch->rho_ccw < branch->rho_cw)
    {
      cw_first = false;
    }
  else
    {
      std::uniform_real_distribution<double> u_dist(0.0,1.0);
      double ur = u_dist(rand_eng);
      if (ur < 0.5) cw_first = false;
    }

  if (cw_first == true)
    {
      while ((((growths[0] < max_growth_cw) && (growths[0] < proposed_r_cw)) ||
	      ((growths[1] < max_growth_ccw) && (growths[1] < proposed_r_ccw))) &&
	     (growths[0] + growths[1] + branch->rho_t < branch->size))
	{
	  if ((growths[0] < max_growth_cw) &&
	      (growths[0] < proposed_r_cw) &&
	      (growths[0] + growths[1] + branch->rho_t < branch->size)) growths[0]++;
	  if ((growths[1] < max_growth_ccw) &&
	      (growths[1] < proposed_r_ccw) &&
	      (growths[0] + growths[1] + branch->rho_t < branch->size)) growths[1]++;
	}
    }
  else
    {
      while ((((growths[0] < max_growth_cw) && (growths[0] < proposed_r_cw)) ||
	      ((growths[1] < max_growth_ccw) && (growths[1] < proposed_r_ccw))) &&
	     (growths[0] + growths[1] + branch->rho_t < branch->size))
	{
	  if ((growths[1] < max_growth_ccw) &&
	      (growths[1] < proposed_r_ccw) &&
	      (growths[0] + growths[1] + branch->rho_t < branch->size)) growths[1]++;
	  if ((growths[0] < max_growth_cw) &&
	      (growths[0] < proposed_r_cw) &&
	      (growths[0] + growths[1] + branch->rho_t < branch->size)) growths[0]++;
	}
    }

  
  return growths;
  
}


// function to perform random growths across active forks
void btree::random_transforms(int r)
{
  int rem = r;
  int potential_size = max_size();
  std::vector<std::string> forks;
  int N_forks, fork_add, fork_rem, r_cw, r_ccw;
  int cycle_rem;
  std::vector<int> partitions;


  // do while units remain  and the total size is the less than the maximum possible size
  while ((rem > 0) &&
	 (total_size() < potential_size))
    {

      // get the set of active forks
      forks = get_active_forks();
      N_forks = count_active_forks();

      partitions.clear();

      cycle_rem = rem;

      std::uniform_int_distribution<int> unif_dist(0,rem);

      // sample partitions
      for (int i=0; i<N_forks; i++)
	{
	  partitions.push_back(unif_dist(rand_eng));
	}

      // sort the partitions
      sort(partitions.begin(),partitions.end());

      // std::cout << "rem = " << rem << std::endl;
      // for (int i=0; i<N_parts; i++)
      // 	{
      // 	  std::cout << partitions[i] << std::endl;
      // 	}

      for (int i=0; i<N_forks; i++)
	{
	  
	  if (partitions[i] > 0)
	    {
	      if (i == 0)
		{
		  fork_add = partitions[i];
		}
	      else if(i == N_forks-1)
		{
		  fork_add = cycle_rem - partitions.back();
		}
	      else
		{
		  fork_add = partitions[i] - partitions[i-1];
		}

	      if (fork_add > 0)
		{
		  std::binomial_distribution<int> binom_dist(fork_add,0.5);
		  r_cw = binom_dist(rand_eng);
		  r_ccw = fork_add - r_cw;
		  fork_rem = grow_at_branch_asym(forks[i],r_cw,r_ccw);
		  rem -= (fork_add - fork_rem);
		}
	    }

	  // std::cout << forks[i] << std::endl;
	  // std::cout << " fork_add = " << fork_add << std::endl;
	  // std::cout << " fork_rem = " << fork_rem << std::endl;

	  // std::cout << " rem = " << rem << std::endl;
	    
	} // end loop over forks
	
    } // end loop over remaining units
  
}


// apply transforms to specified forks using provided rate
void btree::random_transforms_on_forks(std::vector<std::string> forks, double r)
{
  std::poisson_distribution<int> poisson_dist(r);
  int r_cw, r_ccw, r_rem;

  for (size_t i_fork=0; i_fork<forks.size(); i_fork++)
    {
      // sample Poisson distribution for each direction
      r_cw = poisson_dist(rand_eng);
      r_ccw = poisson_dist(rand_eng);
      // execute growth
      // std::cout << forks[i_fork] << std::endl;
      r_rem = grow_at_branch_asym(forks[i_fork],r_cw,r_ccw);
      std::cout << forks[i_fork]
		<< ", rem = "
		<< r_rem
		<< std::endl;
    }
  
}


// function get pointer for a branch given the code
node *btree::get_branch(std::string loc)
{

  int n = loc.length();
  int i;

  node *tar;

  tar = root;

  if (n == 1)
    {
      return tar;
    }
  
  i = 1;
  while ((i < n) && (tar != nullptr))
    {
      tar = parse_dir(tar,loc[i]);
      i += 1;
    }
  return tar;
  
}

node *btree::parse_dir(node *branch, char d)
{
  if (branch->leaf == false)
    {
      if (d == 'l')
	{
	  return branch->left;
	}
      else if (d == 'r')
	{
	  return branch->right;
	}
    }
  return nullptr;
}


// functions to count the total number of leaves
int btree::count_total_leaves()
{
  return count_leaves("m");
}

int btree::count_leaves(std::string loc)
{
  return leaf_counter(get_branch(loc));
}

int btree::leaf_counter(node *branch)
{
  int leaf_count = 0;
  if (branch->leaf == true)
    {
      leaf_count += 1;
    }
  else
    {
      leaf_count += leaf_counter(branch->left);
      leaf_count += leaf_counter(branch->right);
    }
  return leaf_count;
}


// functions to calculate the total number of forks
int btree::count_total_forks()
{
  return count_total_forks("m");
}

int btree::count_total_forks(std::string loc)
{
  return total_fork_counter(get_branch(loc));
}

int btree::total_fork_counter(node *branch)
{
  int fork_count = 0;
  if (branch->leaf == false)
    {
      fork_count += 1;
      fork_count += total_fork_counter(branch->left);
      fork_count += total_fork_counter(branch->right);
    }
  return fork_count;
}


// counting functions for completed forks
int btree::count_completed_forks()
{
  return count_completed_forks("m");
}

int btree::count_completed_forks(std::string loc)
{
  return completed_fork_counter(get_branch(loc));
}

int btree::completed_fork_counter(node *branch)
{
  int fork_count = 0;
  if (branch->leaf == false)
    {
      if (branch->complete == true) fork_count += 1;
      fork_count += completed_fork_counter(branch->left);
      fork_count += completed_fork_counter(branch->right);
    }
  return fork_count;
}


// counting functions for active forks
int btree::count_active_forks()
{
  return count_active_forks("m");
}

int btree::count_active_forks(std::string loc)
{
  return active_fork_counter(get_branch(loc));
}

int btree::active_fork_counter(node *branch)
{
  int fork_count = 0;
  if (branch->leaf == false)
    {
      if (branch->complete == false) fork_count += 1;
      fork_count += active_fork_counter(branch->left);
      fork_count += active_fork_counter(branch->right);
    }
  return fork_count;
}


// functions to get labels for leaves
std::vector<std::string> btree::get_leaves()
{
  std::vector<std::string> leaves;
  if (count_total_leaves() > 0)
    {
      leaves.push_back("m");
      leaves = traverse_leaves(leaves,get_branch("m"));
    }
  return leaves;
}

std::vector<std::string> btree::traverse_leaves(std::vector<std::string> leaves, node *branch)
{
  if (branch->leaf == false)
    {
      std::string temp_parent = leaves.back();
      leaves.back() += "l";
      leaves = traverse_leaves(leaves,branch->left);
      leaves.push_back(temp_parent+"r");
      leaves = traverse_leaves(leaves,branch->right);
    }
  return leaves;
}


// functions to get labels for completed forks
std::vector<std::string> btree::get_completed_forks()
{
  std::vector<std::string> forks;
  if (count_completed_forks() > 0)
    {
      forks = traverse_completed_forks(forks,0,get_branch("m"));
    }
  return forks;
}

std::vector<std::string> btree::traverse_completed_forks(std::vector<std::string> forks, int i, node *branch)
{
  
  // case for root node
  if (forks.size() == 0)
    {
      if ((branch->leaf == false) &&
	  (branch->complete == true))
	{ 
	  forks.push_back("m");
	}
    }

  // case for other nodes
  if ((branch->leaf == false) &&
      (branch->complete == true))
    {
      forks.push_back(forks[i]+"l");
      forks = traverse_completed_forks(forks,forks.size()-1,branch->left);
      forks.push_back(forks[i]+"r");
      forks = traverse_completed_forks(forks,forks.size()-1,branch->right);
    }
  else
    {
      forks.pop_back();
    }
  return forks;
}


// functions to get labels for active forks
std::vector<std::string> btree::get_active_forks()
{
  std::vector<std::string> forks;
  if (count_active_forks() > 0)
    {
      forks = traverse_active_forks(forks,"m",get_branch("m"));
    }
  return forks;
}

std::vector<std::string> btree::traverse_active_forks(std::vector<std::string> forks, std::string f, node *branch)
{
  
  if (branch->leaf == false)
    { 
      if (branch->complete == false)
	{ 
	  forks.push_back(f);
	}
      forks = traverse_active_forks(forks,f+"l",branch->left);
      forks = traverse_active_forks(forks,f+"r",branch->right);
    }
  return forks;
}


// calculate the total size of the system
int btree::total_size()
{
  return branch_size(get_branch("m"));
}

// calculate the scaled size of the system
double btree::scaled_size()
{
  double s = 1.0*total_size();
  return s/single_size();
}


// calculate size of branch
int btree::branch_size(node *branch)
{
  int s = 0;

  if (branch->parent == nullptr)
    {
      s += branch->size;
    }
  
  if (branch->leaf == false)
    {
      s += branch->rho_t;
      s += branch_size(branch->left);
      s += branch_size(branch->right);
    }
  
  return s;
}


// calculate size of a separable branch
int btree::separable_branch_size(node *branch)
{
  int s = 0;

  if (branch->parent == nullptr)
    {
      s += branch->size;
    }
  else
    {
      s += branch->parent->rho_t;
      if (branch->leaf == false)
	{
	  s -= branch->rho_t;
	  s += separable_branch_size(branch->left);
	  s += separable_branch_size(branch->right);
	}
    }
  
  return s;
}


// calculate the maximum possible size of the system
int btree::max_size()
{
  return (root->size)*count_total_leaves();
}

// get the size of a single genome
int btree::single_size()
{
  return root->size;
}


// function to partition units among branches
void btree::solve_topology()
{
  int temp_start, target_size;
  std::string free_leaf;
  std::string query_leaf;
  node *free_branch; // free branch whose topology is being determined
  node *query_branch; // branch queried for sizing of free branch
  node *link_branch; // branch that free branch is linked to

  int N_leaves = count_total_leaves();
  std::vector<std::string> leaves = get_leaves();

  bool match;
  int j_q, j_f;

  int mid;
  //int mid, d_start, d_end;

  temp_start = 0;

  for (int i=0; i<N_leaves; i++)
    {

      free_leaf = leaves[i];
      free_branch = get_branch(free_leaf);

      // first free branch is always circular
      if (i == 0)
	{
	  
	  target_size = free_branch->size;
	  
	}
      // find query branch supporting free branch
      else
	{

	  match = false;
	  query_leaf = leaves[i-1];

	  j_f = free_leaf.length() - 1;
	  while ((match == false) && (j_f >= 0))
	    {

	      j_q = query_leaf.length() - 1;
	      while ((match == false) && (j_q >= 0))
		{

		  if (free_leaf.substr(0,j_f) == query_leaf.substr(0,j_q))
		    {
		      query_branch = get_branch(query_leaf.substr(0,j_q));
		      match = true;
		    }

		  j_q -= 1;
	      
		}

	      j_f -= 1;
	      
	    }

	  // set the target size from query branch
	  target_size = query_branch->rho_t;

	  // set the link_branch to the query branch
	  link_branch = query_branch;

	  // descend to leftmost leaf of link branch
	  while (link_branch->leaf == false)
	    {
	      link_branch = link_branch->left;
	    }
	  
	}

      // start constructing the topology


      if (target_size > 0)
	{

	  // set the start, end, and mid
	  free_branch->topo.start = temp_start;
	  free_branch->topo.end = temp_start + target_size - 1;
	  // free_branch->topo.mid = (free_branch->topo.start + free_branch->topo.end)/2;

	  // if on the first branch or the query branch is complete, then set midpoint using size
	  if ((i == 0) ||
	      (query_branch->complete == true)) // query branch will be undefined on first iteration
	    {
	      // set midpoint
	      free_branch->topo.mid = free_branch->topo.start + target_size/2;
	      free_branch->topo.start_link = free_branch->topo.end;
	      // create a circle
	      free_branch->topo.end_link = free_branch->topo.start;
	    }
	  // otherwise determine the midpoint and create a theta structure
	  else
	    {

	      // align midpoint with link branch
	      mid = link_branch->topo.mid;
	      
	      if (query_branch->rho_cw == 0)
		{
		  // determine the midpoint
		  free_branch->topo.mid = free_branch->topo.end;
		  // determine the start link
		  if ((mid - query_branch->rho_ccw) > link_branch->topo.start)
		    {
		      free_branch->topo.start_link = mid - query_branch->rho_ccw;
		    }
		  else
		    {
		      free_branch->topo.start_link = link_branch->topo.end +
			((mid - query_branch->rho_ccw) - link_branch->topo.start); // BRG - caution
		    }
		  // determine the end link
		  free_branch->topo.end_link = mid + 1;
		}
	      else if (query_branch->rho_ccw == 0)
		{
		  // determine the midpoint
		  free_branch->topo.mid = free_branch->topo.start;
		  // determine the start link
		  free_branch->topo.start_link = mid - 1;
		  // determine the end link
		  if ((mid + query_branch->rho_cw) <= link_branch->topo.end)
		    {
		      free_branch->topo.end_link = mid + query_branch->rho_cw;
		    }
		  else
		    {
		      free_branch->topo.end_link = (link_branch->topo.start) +
			((mid + query_branch->rho_cw) - link_branch->topo.end) - 1; // BRG - caution
		    }
		}
	      else
		{
		  // determine the midpoint
		  free_branch->topo.mid = free_branch->topo.start + query_branch->rho_ccw;

		  // determine the start link
		  if ((mid - query_branch->rho_ccw - 1) > link_branch->topo.start)
		    {
		      free_branch->topo.start_link = mid - query_branch->rho_ccw - 1;
		    }
		  else
		    {
		      // free_branch->topo.start_link = link_branch->topo.end +
		      // 	((mid - query_branch->rho_ccw) - link_branch->topo.start); // BRG - caution

		      free_branch->topo.start_link = link_branch->topo.start;
		    }

		  // determine the end link
		  if ((mid + query_branch->rho_cw) <= link_branch->topo.end)
		    {
		      free_branch->topo.end_link = mid + query_branch->rho_cw;
		    }
		  else
		    {
		      free_branch->topo.end_link = (link_branch->topo.start) +
			((mid + query_branch->rho_cw) - link_branch->topo.end) - 1; // BRG - caution
		    }
		}

	      if (target_size == free_branch->size - 1) free_branch->topo.end_link = free_branch->topo.start_link;
	    }

	} // end if for target_size

      	  // increase the starting location by the added size
      temp_start += target_size;
      
    } // end loop over free branches
  
}


// function used to dump the topology to a file
void btree::dump_topology(std::string topo_filename, int idx)
{
  std::fstream topo_file;
  
  node *topo_branch;

  topo_file.open(topo_filename, std::ios::out);

  if (!topo_file.is_open())
    {
      std::cout << "ERROR: file not opened in dump_topology" << std::endl;
    }
  else
    {

      topo_file << "size=" << total_size() <<  std::endl;      
  
      for (std::string leaf: get_leaves())
	{
	  
	  topo_branch = get_branch(leaf);

	  topo_file << leaf
		    << "("
		    << (topo_branch->topo.end-topo_branch->topo.start+1)
		    << ")"
		    << "," << (topo_branch->topo.start_link + idx)
		    << "," << (topo_branch->topo.start + idx)
		    << "," << (topo_branch->topo.mid + idx)
		    << "," << (topo_branch->topo.end + idx)
		    << "," << (topo_branch->topo.end_link + idx)
		    << std::endl;
      
	}

      topo_file.close();
      
    }

}


// get the topology of a single leaf
theta_topo btree::get_leaf_topo(std::string loc)
{
  node *branch;

  branch = get_branch(loc);

  return branch->topo;
}


// function used to dump the fork partitions to a file
void btree::dump_fork_partitions(std::string fork_partitions_filename, int idx)
{

  std::vector<fork_partition> f_ps = get_all_fork_partitions();
  
  std::fstream f_ps_file;
  std::string temp_line;
  int N_total;
  
  f_ps_file.open(fork_partitions_filename, std::ios::out);

  if (!f_ps_file.is_open())
    {
      std::cout << "ERROR: file not opened in dump_fork_partitions" << std::endl;
    }
  else
    {

      // write the total number of forks that were partitioned about
      f_ps_file << "N_forks=" << f_ps.size() << std::endl;
      // write the indexing convention
      f_ps_file << "idx=" << idx << std::endl;

      for (fork_partition f_p : f_ps)
	{

	  temp_line = f_p.fork + ",";
	  temp_line += std::to_string(f_p.left_monos.size()) + ",";

	  N_total = 0;
	  for (size_t i=0; i<f_p.left_monos.size(); i++)
	    N_total += f_p.left_monos[i].N;
	  temp_line += std::to_string(N_total) + ",";
	  
	  temp_line += std::to_string(f_p.right_monos.size()) + ",";

	  N_total = 0;
	  for (size_t i=0; i<f_p.right_monos.size(); i++)
	    N_total += f_p.right_monos[i].N;
	  temp_line += std::to_string(N_total);
	  
	  f_ps_file << temp_line << std::endl;

	  for (size_t i=0; i<f_p.left_monos.size(); i++)
	    {

	      if (f_p.left_monos[i].wrapped == false)
		{
		  temp_line = "\tnw,";
		}
	      else
		{
		  temp_line = "\tw,";
		}

	      temp_line += std::to_string(f_p.left_monos[i].N) + ",";
	      temp_line += std::to_string(f_p.left_monos[i].ll) + ",";
	      temp_line += std::to_string(f_p.left_monos[i].mid_ll) + ",";
	      temp_line += std::to_string(f_p.left_monos[i].mid_ul) + ",";
	      temp_line += std::to_string(f_p.left_monos[i].ul);

	      f_ps_file << temp_line << std::endl;
	      
	    }

	  for (size_t i=0; i<f_p.right_monos.size(); i++)
	    {

	      if (f_p.right_monos[i].wrapped == false)
		{
		  temp_line = "\tnw,";
		}
	      else
		{
		  temp_line = "\tw,";
		}

	      temp_line += std::to_string(f_p.right_monos[i].N) + ",";
	      temp_line += std::to_string(f_p.right_monos[i].ll) + ",";
	      temp_line += std::to_string(f_p.right_monos[i].mid_ll) + ",";
	      temp_line += std::to_string(f_p.right_monos[i].mid_ul) + ",";
	      temp_line += std::to_string(f_p.right_monos[i].ul);

	      f_ps_file << temp_line << std::endl;
	      
	    }
	}

      f_ps_file.close();
      
    }
}


// get all fork partitions
std::vector<fork_partition> btree::get_all_fork_partitions()
{
  std::vector<fork_partition> f_ps;
  fork_partition f_p;

  // get the total set of forks (completed and active)
  
  std::vector<std::string> total_forks;
  
  for (std::string fork : get_completed_forks())
    {
      total_forks.push_back(fork);
    }
  for (std::string fork : get_active_forks())
    {
      total_forks.push_back(fork);
    }


  // iterate over the forks and get the partition for each
  for (std::string fork : total_forks)
    {
      f_p = get_fork_partition(fork);
      f_ps.push_back(f_p);
    }

  return f_ps;
  
}


// get the partition about a fork
fork_partition btree::get_fork_partition(std::string loc)
{

  fork_partition f_p;

  node *forked_branch = get_branch(loc);
  f_p.fork = loc;

  // determine the leaves belonging to the left and right branches
  std::string l_loc, r_loc;
  l_loc = loc + 'l';
  r_loc = loc + 'r';
  std::vector<std::string> l_leaves, r_leaves;
  l_leaves.push_back(l_loc);
  l_leaves = traverse_leaves(l_leaves,get_branch(l_loc));
  r_leaves.push_back(r_loc);
  r_leaves = traverse_leaves(r_leaves,get_branch(r_loc));

  // prepare the monomer ranges
  theta_topo temp_topo;
  mono_range temp_m_r;
  std::string leaf;

  // add the monomer ranges for the left leaves
  for (size_t i_leaf=0; i_leaf<l_leaves.size(); i_leaf++)
    {
      
      // initialize the mono_range
      temp_m_r.wrapped = false;
      temp_m_r.mid_ll = -1;
      temp_m_r.mid_ul = -1;

      // get the topology of the leaf
      leaf = l_leaves[i_leaf];
      temp_topo = get_leaf_topo(leaf);

      // add the monomer range
      temp_m_r.ll = temp_topo.start;
      temp_m_r.ul = temp_topo.end;
      temp_m_r.N = temp_m_r.ul - temp_m_r.ll + 1;
      f_p.left_monos.push_back(temp_m_r);

    }

  // add the monomer ranges for the right leaves
  for (size_t i_leaf=0; i_leaf<r_leaves.size(); i_leaf++)
    {
      
      // initialize the mono_range
      temp_m_r.wrapped = false;
      temp_m_r.mid_ll = -1;
      temp_m_r.mid_ul = -1;

      // get the topology of the leaf
      leaf = r_leaves[i_leaf];
      temp_topo = get_leaf_topo(leaf);

      // add the monomer range
      temp_m_r.ll = temp_topo.start;
      temp_m_r.ul = temp_topo.end;
      temp_m_r.N = temp_m_r.ul - temp_m_r.ll + 1;
      f_p.right_monos.push_back(temp_m_r);

    }

  int N;

  // perform the correction for the overcounting on the left side
  if (forked_branch->complete == false)
    {
      
      temp_topo = get_leaf_topo(r_leaves[0]);

      // replication has proceeded past Ter but forks have not met
      if (temp_topo.start_link >= temp_topo.end_link)
	{
	  
	  f_p.left_monos[0].wrapped = true;
	  // extended past Ter in ccw direction
	  if (forked_branch->rho_cw < forked_branch->rho_ccw)
	    {
	      f_p.left_monos[0].mid_ll = temp_topo.end_link - 1;
	      if (temp_topo.start_link < f_p.left_monos[0].ul)
		{
		  f_p.left_monos[0].mid_ul = temp_topo.start_link + 1;
		}
	      else
		{
		  f_p.left_monos[0].ul = -2;
		  f_p.left_monos[0].mid_ul = -1;
		}
	    }
	  // extended past Ter in cw direction
	  else
	    {
	      f_p.left_monos[0].mid_ul = temp_topo.start_link + 1;
	      if (temp_topo.end_link > f_p.left_monos[0].ll)
		{
		  f_p.left_monos[0].mid_ll = temp_topo.end_link - 1;
		}
	      else
		{
		  f_p.left_monos[0].ll = -1;
		  f_p.left_monos[0].mid_ll = -2;
		}
	    }

	  // count the number of monomers in the new ranges
	  N = 0;
	  for (int i=f_p.left_monos[0].ll; i<(f_p.left_monos[0].mid_ll+1); i++)
	    {
	      N += 1;
	    }
	  for (int i=f_p.left_monos[0].mid_ul; i<(f_p.left_monos[0].ul+1); i++)
	    {
	      N += 1;
	    }
	  f_p.left_monos[0].N = N;
	  
	}
      // replication has not proceeded past the Ter
      else
	{
	  f_p.left_monos[0].ll = temp_topo.start_link + 1;
	  f_p.left_monos[0].ul = temp_topo.end_link - 1;
	  f_p.left_monos[0].N = f_p.left_monos[0].ul - f_p.left_monos[0].ll + 1;
	}

    }

  return f_p;
}


// function to read chromosome regions from file
std::vector<chromo_region> btree::read_regions(std::string rg_filename, int idx)
{
  std::fstream rg_file;

  chromo_region c_r;
  std::vector<chromo_region> c_rs;

  std::string line;
  
  std::string rg_delim = ",";
  int delim;
  
  rg_file.open(rg_filename, std::ios::in);

  if (!rg_file.is_open())
    {
      std::cout << "ERROR: file not opened in read_regions" << std::endl;
    }
  else
    {
      while (1)
	{
	  rg_file >> line;
	  if (rg_file.eof()) break;

	  // region name
	  delim = line.find(rg_delim);
	  c_r.name = line.substr(0,delim);
	  line.erase(0,delim+1);

	  // start of region (inclusive)
	  delim = line.find(rg_delim);
	  c_r.start = stoi(line.substr(0,delim)) - idx;
	  line.erase(0,delim+1);

	  // end of region (inclusive)
	  c_r.end = stoi(line) - idx;

	  // initialize count to zero
	  c_r.count = 0;

	  c_rs.push_back(c_r);
	  
	}

      rg_file.close();
      
    }

  return c_rs;
  
}


// function to update region counts
void btree::update_region_counts(std::vector<chromo_region> &c_rs)
{

  int offset;
  int circ_size, exists_size;
  int start_circ, end_circ;
  node *count_branch;

  // zero counts before summation over leaves
  for (chromo_region &c_r: c_rs)
    {
      c_r.count = 0;
    }

  // iterate over leaves
  for (std::string leaf: get_leaves())
    {
	  
      count_branch = get_branch(leaf);

      // calculate appropriate offset from topology
      circ_size = count_branch->size;
      exists_size = count_branch->topo.end - count_branch->topo.start + 1;
      offset = count_branch->topo.mid - count_branch->topo.start + 1;
      // exists_size = count_branch->parent->rho_t;
      // offset = exists_size/2;
      // offset = count_branch->parent->rho_ccw; 

      // std::cout << exists_size << std::endl;
      // std::cout << offset << std::endl;

      // iterate over regions
      for (chromo_region &c_r: c_rs)
	{

	  start_circ = (c_r.start + offset)%circ_size;
	  end_circ = (c_r.end + offset)%circ_size;

	  // std::cout << c_r.name << " " << start_circ << " " << end_circ << std::endl;

	  if ((start_circ >= 0) &&
	      (end_circ >= 0) &&
	      (start_circ < exists_size) &&
	      (end_circ < exists_size))
	    {

	      c_r.count += 1;

	    }

	}

    }
  
}


// function used to dump the regions and counts to a file
void btree::dump_regions(std::string rg_filename, std::vector<chromo_region> c_rs)
{
  std::fstream rg_file;

  rg_file.open(rg_filename, std::ios::out);

  if (!rg_file.is_open())
    {
      std::cout << "ERROR: file not opened in dump_regions" << std::endl;
    }
  else
    {
  
      for (chromo_region c_r: c_rs)
	{

	  rg_file << c_r.name
		  << "," << c_r.start
		  << "," << c_r.end
		  << ":" << c_r.count
		  << std::endl;
      
	}

      rg_file.close();

    }
  
}

void btree::dump_CG_map(std::string CG_filename, int idx, CG_map &m)
{
  std::fstream CG_file;

  CG_file.open(CG_filename, std::ios::out);

  if (!CG_file.is_open())
    {
      std::cout << "ERROR: file not opened in dump_CG_map" << std::endl;
    }
  else
    {

      CG_file << "N_base = "
	      << m.N_base
	      << "\nN = "
	      << m.N
	      << "\nf_CG = "
	      << m.f_CG
	      << "\nN_base_CG = "
	      << m.N_base_CG
	      << "\nN_CG = "
	      << m.N_CG
	      << std::endl;

      CG_file << "\n\nN_leaves = "
	      << m.CG_leaves.size()
	      << std::endl;

      for (CG_leaf leaf: m.CG_leaves)
	{

	  CG_file << leaf.leaf
		  << ","
		  << leaf.start + idx
		  << ","
		  << leaf.end + idx
		  << std::endl;
      
	}

      CG_file << "\nID, base-ID, min, max" << std::endl;
      
  
      for (CG_locus l: m.loci)
	{

	  CG_file << l.CG + idx
		  << ","
		  << l.bCG + idx
		  << ","
		  << l.start + idx
		  << ","
		  << l.end + idx
		  << std::endl;
      
	}

      CG_file.close();

    }

}


// function to create a coarse-graining for the entire system
CG_map btree::update_CG_map(int f_CG)
{

  CG_map m;
  CG_leaf temp_CG_leaf;

  int temp_lb, temp_ub, temp_mid, mid_offset;

  std::vector<std::string> leaves = get_leaves();

  m.f_CG = f_CG;
  m.N_base = root->size;
  m.N = total_size();

  for (size_t i_leaf=0; i_leaf<leaves.size(); i_leaf++)
    {

      // get the current number of loci
      temp_lb = m.loci.size();

      // add loci based on a centered CG map
      centered_CG_map(temp_mid,m.loci,get_branch(leaves[i_leaf]),f_CG);

      // get the updated number of loci
      temp_ub = m.loci.size();

      temp_CG_leaf.leaf = leaves[i_leaf];
      temp_CG_leaf.start = temp_lb;
      temp_CG_leaf.end = temp_ub - 1;
      temp_CG_leaf.mid = temp_mid + temp_lb;

      if (i_leaf == 0)
	{
	  // set the base number of CG to the first leaf's number of loci
	  m.N_base_CG = m.loci.size();
	}
      else
	{

	  // calculate the mid offset to correctly align Oris
	  mid_offset = (m.CG_leaves[0].mid - m.CG_leaves[0].start) - 
	    (temp_CG_leaf.mid - temp_CG_leaf.start);

	  // correct loci indices for subsequent leaves
	  for (int j=temp_CG_leaf.start; j<temp_CG_leaf.end+1; j++)
	    {
	      m.loci[j].CG += (m.loci[temp_CG_leaf.start-1].CG + 1);
	      m.loci[j].bCG += mid_offset;
	      m.loci[j].start += (m.loci[temp_CG_leaf.start-1].end + 1);
	      m.loci[j].end += (m.loci[temp_CG_leaf.start-1].end + 1);
	    }
	  
	}
      
      m.CG_leaves.push_back(temp_CG_leaf);
      
    }

  m.N_CG = m.loci.size();

  return m;
}


// function to create a centered coarse-graining for a single branch
void btree::centered_CG_map(int &mid, std::vector<CG_locus> &loci, node *branch, int f_CG)
{

  CG_locus l;
  // int N = branch->topo.end - branch->topo.start + 1;
  // int ori_idx_m = branch->topo.mid - branch->topo.start;
  // int ori_idx_m = branch->topo.mid - branch->topo.start - 1;
  // int ori_idx_p = ori_idx_m + 1;
  int ori_idx_m = branch->topo.mid - branch->topo.start;
  int ori_idx_p = branch->topo.end - branch->topo.mid + 1;

  int N_cum, N_CG_cum;

  int CG_rem_m, CG_rem_p, N_CG_m, N_CG_p, N_CG;

  CG_rem_m = ori_idx_m%f_CG;
  CG_rem_p = ori_idx_p%f_CG;

  // N_CG_m = (ori_idx_m+1)/f_CG;
  N_CG_m = ori_idx_m/f_CG;
  if (CG_rem_m > 0) N_CG_m++;
  // N_CG_p = (N - ori_idx_p)/f_CG;
  N_CG_p = ori_idx_p/f_CG;
  if (CG_rem_p > 0) N_CG_p++;

  N_CG = N_CG_m + N_CG_p;
  mid = N_CG_m + 1;

  N_cum = 0;
  N_CG_cum = 0;

  while (N_CG_cum < N_CG)
    {

      // conditional for remainder beads at start
      if ((N_CG_cum == 0) &&
	  (CG_rem_m > 0))
	{

	  l.CG = N_CG_cum;
	  l.bCG = N_CG_cum;
	  l.start = N_cum;
	  l.end = N_cum + CG_rem_m - 1;

	  N_CG_cum += 1;
	  N_cum += CG_rem_m;

	  loci.push_back(l);

	  continue;
	  
	}

      // conditional for remainder beads at end
      if ((N_CG_cum == N_CG - 1) &&
	  (CG_rem_p > 0))
	{

	  l.CG = N_CG_cum;
	  l.bCG = N_CG_cum;
	  l.start = N_cum;
	  l.end = N_cum + CG_rem_p - 1;

	  N_CG_cum += 1;
	  N_cum += CG_rem_p;

	  loci.push_back(l);

	  continue;
	  
	}

      l.CG = N_CG_cum;
      l.bCG = N_CG_cum;
      l.start = N_cum;
      l.end = N_cum + f_CG - 1;

      N_CG_cum += 1;
      N_cum += f_CG;

      loci.push_back(l);
      

    }
  
}


// function to prepare the bonds
void btree::prepare_bonds(int **&c, int *&t, int &N, int idx)
{

  node *topo_leaf;
  
  N = total_size() + count_active_forks();

  t = new int[N];

  c = new int*[N];
  for (int i=0; i<N; i++)
    {
      c[i] = new int[2];
    }

  int i_bond = 0;
  // iterate over leaves in system
  std::vector<std::string>leaves = get_leaves();
  for (int i_leaf=0; i_leaf<count_total_leaves(); i_leaf++)
    {
      
      topo_leaf = get_branch(leaves[i_leaf]);

      // test if leaf is complete
      if ((topo_leaf->topo.start_link == topo_leaf->topo.end) &&
	  (topo_leaf->topo.end_link == topo_leaf->topo.start))
	{
	  for (int i=0; i<topo_leaf->size; i++)
	    {
	      t[i_bond] = 1;
	      c[i_bond][0] = topo_leaf->topo.start + i;
	      c[i_bond][1] = topo_leaf->topo.start + (i+1)%(topo_leaf->size);
	      i_bond += 1;
	    }
	}
      else
	{
	  // bond at start
	  t[i_bond] = 1;
	  c[i_bond][0] = topo_leaf->topo.start_link;
	  c[i_bond][1] = topo_leaf->topo.start;
	  i_bond += 1;

	  // bonds in middle
	  for (int i=topo_leaf->topo.start; i<topo_leaf->topo.end; i++)
	    {
	      t[i_bond] = 1;
	      c[i_bond][0] = i;
	      c[i_bond][1] = i + 1;
	      i_bond += 1;
	    }

	  // bond at end
	  t[i_bond] = 1;
	  c[i_bond][0] = topo_leaf->topo.end;
	  c[i_bond][1] = topo_leaf->topo.end_link;
	  i_bond += 1;
	}
      
    }
  
  // apply indexing convention
  for (int i=0; i<N; i++)
    {
      for (int j=0; j<2; j++)
	{
	  c[i][j] += idx;
	}
    }
}


// function to prepare the angles
void btree::prepare_angles(int **&c, int *&t, int &N, int idx)
{
  
  node *topo_leaf = nullptr;
  node *supp_leaf = nullptr;
  
  N = 2*(total_size() + count_active_forks());

  t = new int[N];

  c = new int*[N];
  for (int i=0; i<N; i++)
    {
      c[i] = new int[3];
    }

  int i_angle = 0;
  // iterate over leaves in system
  std::vector<std::string>leaves = get_leaves();
  for (int i_leaf=0; i_leaf<count_total_leaves(); i_leaf++)
    {
      
      topo_leaf = get_branch(leaves[i_leaf]);

      // test if leaf is complete
      if ((topo_leaf->topo.start_link == topo_leaf->topo.end) &&
	  (topo_leaf->topo.end_link == topo_leaf->topo.start))
	{
	  for (int i=0; i<topo_leaf->size; i++)
	    {
	      t[i_angle] = 1;
	      c[i_angle][0] = topo_leaf->topo.start + (i+topo_leaf->size-1)%(topo_leaf->size);
	      c[i_angle][1] = topo_leaf->topo.start + i;
	      c[i_angle][2] = topo_leaf->topo.start + (i+1)%(topo_leaf->size);
	      i_angle += 1;

	      t[i_angle] = 2;
	      for (int j=0; j<3; j++)
		{
		  c[i_angle][j] = c[i_angle-1][j];
		}
	      i_angle += 1;
	    }
	}
      else
	{

	  // determine chromosome containing start_link and end_link, thereby supporting the current topo_leaf
	  for (int j_leaf=0; j_leaf<i_leaf; j_leaf++)
	    {
	      supp_leaf = get_branch(leaves[j_leaf]);

	      if ((topo_leaf->topo.start_link >= supp_leaf->topo.start) &&
		  (topo_leaf->topo.end_link <= supp_leaf->topo.end)) break;
	    }

	  // correct angles centered about start_link on supp_leaf
	  for (int j_angle=0; j_angle<i_angle; j_angle++)
	    {
	      if (c[j_angle][1] == topo_leaf->topo.start_link)
		{
		  if (t[j_angle] == 1)
		    {
		      t[j_angle] = 3;
		    }
		  else if(t[j_angle] == 2)
		    {
		      t[j_angle] = 4;
		    }
		}
	    }

	  // angle centered about start_link on topo_leaf
	  t[i_angle] = 3;
	  if (topo_leaf->topo.start_link == supp_leaf->topo.start)
	    {
	      c[i_angle][0] = supp_leaf->topo.start_link;
	    }
	  else
	    {
	      c[i_angle][0] = topo_leaf->topo.start_link - 1;
	    }
	  c[i_angle][1] = topo_leaf->topo.start_link;
	  c[i_angle][2] = topo_leaf->topo.start;
	  i_angle += 1;

	  t[i_angle] = 4;
	  for (int j=0; j<3; j++)
	    {
	      c[i_angle][j] = c[i_angle-1][j];
	    }
	  i_angle += 1;

	  // angle at start
	  t[i_angle] = 1;
	  c[i_angle][0] = topo_leaf->topo.start_link;
	  c[i_angle][1] = topo_leaf->topo.start;
	  c[i_angle][2] = std::min(topo_leaf->topo.start+1,topo_leaf->topo.end);
	  i_angle += 1;

	  t[i_angle] = 2;
	  for (int j=0; j<3; j++)
	    {
	      c[i_angle][j] = c[i_angle-1][j];
	    }
	  i_angle += 1;

	  // angles in middle
	  for (int i=topo_leaf->topo.start+1; i<topo_leaf->topo.end-1; i++)
	    {
	      t[i_angle] = 1;
	      c[i_angle][0] = i - 1;
	      c[i_angle][1] = i;
	      c[i_angle][2] = i + 1;
	      i_angle += 1;

	      t[i_angle] = 2;
	      for (int j=0; j<3; j++)
		{
		  c[i_angle][j] = c[i_angle-1][j];
		}
	      i_angle += 1;
	    }

	  // angle at end
	  t[i_angle] = 1;
	  c[i_angle][0] = std::max(topo_leaf->topo.end-1,topo_leaf->topo.start);
	  c[i_angle][1] = topo_leaf->topo.end;
	  c[i_angle][2] = topo_leaf->topo.end_link;
	  i_angle += 1;
	  
	  t[i_angle] = 2;
	  for (int j=0; j<3; j++)
	    {
	      c[i_angle][j] = c[i_angle-1][j];
	    }
	  i_angle += 1;

	  // correct angles centered about end_link on supp_leaf
	  for (int j_angle=0; j_angle<i_angle; j_angle++)
	    {
	      if (c[j_angle][1] == topo_leaf->topo.end_link)
		{
		  if (t[j_angle] == 1)
		    {
		      t[j_angle] = 3;
		    }
		  else if(t[j_angle] == 2)
		    {
		      t[j_angle] = 4;
		    }
		}
	    }

	  // angle centered about end_link on topo_leaf
	  t[i_angle] = 3;
	  c[i_angle][0] = topo_leaf->topo.end;
	  c[i_angle][1] = topo_leaf->topo.end_link;
	  if (topo_leaf->topo.end_link == supp_leaf->topo.end)
	    {
	      c[i_angle][2] = supp_leaf->topo.end_link;
	    }
	  else
	    {
	      c[i_angle][2] = topo_leaf->topo.end_link - 1;
	    }
	  i_angle += 1;

	  t[i_angle] = 4;
	  for (int j=0; j<3; j++)
	    {
	      c[i_angle][j] = c[i_angle-1][j];
	    }
	  i_angle += 1;
	  
	}
    }

  // apply indexing convention
  for (int i=0; i<N; i++)
    {
      for (int j=0; j<3; j++)
	{
	  c[i][j] += idx;
	}
    }  
}


// prepare the types
void btree::prepare_types(int *&t, int &N, int base_type)
{
  node *topo_leaf;

  t = new int[N];

  // iterate over leaves in system
  std::vector<std::string>leaves = get_leaves();
  for (int i_leaf=0; i_leaf<count_total_leaves(); i_leaf++)
    {
      
      topo_leaf = get_branch(leaves[i_leaf]);

      // unmodified monomers
      for (int i=topo_leaf->topo.start; i<topo_leaf->topo.end+1; i++)
	{
	  t[i] = base_type;
	}

      // add oris at the midpoints
      t[topo_leaf->topo.mid] = base_type + 1;

      // create a ter if the leaf is complete
      if ((topo_leaf->topo.start_link == topo_leaf->topo.end) &&
	  (topo_leaf->topo.end_link == topo_leaf->topo.start))
	{
	  t[topo_leaf->topo.start] = base_type + 2;
	}
      else // create forks
	{
	  t[topo_leaf->topo.start_link] = base_type + 3;
	  t[topo_leaf->topo.end_link] = base_type + 3;
	}

    }
}


// function used by destructor
void btree::destroy_tree(node *branch)
{
  if (branch->leaf == false)
    {
      destroy_tree(branch->left);
      destroy_tree(branch->right);
    }
  delete branch;
}

void btree::destroy_tree()
{
  if (root != nullptr)
    {
      destroy_tree(root);
      reset_root();
    }
}


// function to print branches of the tree
void btree::print_branch(node *branch)
{

  std::string gen_offset;

  gen_offset = "";
  for (int i=0; i<branch->gen; i++)
    {
      gen_offset += "  ";
    }
  gen_offset += "| ";
  
  std::cout << gen_offset
	    << "generation = "
	    << branch->gen
	    << std::endl;
  
  if (branch->leaf == false)
    {
      int max_size_cw, max_size_ccw, max_size_t;
      if (branch->parent != nullptr)
	{
	  max_size_t = branch->parent->rho_t;
	  max_size_cw = branch->parent->rho_cw;
	  max_size_ccw = branch->parent->rho_ccw;
	}
      else
	{
	  max_size_t = branch->size;
	  // max_size_cw = max_size_t - branch->rho_ccw;
	  // max_size_ccw = max_size_t - branch->rho_cw;
	  max_size_cw = max_size_t;
	  max_size_ccw = max_size_t;
	}
      max_size_cw -= branch->rho_ccw;
      max_size_ccw -= branch->rho_cw;
      
      std::cout << gen_offset
	   << "rho_t = "
	   << branch->rho_t
	   << "/" << max_size_t
	   << ", rho_cw = "
	   << branch->rho_cw
	   << "/" << max_size_cw
	   << ", rho_ccw = "
	   << branch->rho_ccw
	   << "/" << max_size_ccw
	   << std::endl;
    }


  std::cout << gen_offset
       << "start = "
       << branch->topo.start
       << ", mid = "
       << branch->topo.mid
       << ", end = "
       << branch->topo.end
       << std::endl;
  std::cout << gen_offset
       << "start_link = "
       << branch->topo.start_link
       << ", end_link = "
       << branch->topo.end_link
       << std::endl;

  if (branch->leaf == false)
    {
      std::cout << gen_offset
		<< "left branch"
		<< std::endl;
      print_branch(branch->left);
      std::cout << gen_offset
		<< "right branch"
		<< std::endl;
      print_branch(branch->right);
    }

}


// function to print the entire tree
void btree::print_tree()
{
  
  if (root != nullptr)
    {

      int N_leaves = count_total_leaves();
      int N_forks = count_total_forks();
      int N_completed_forks = count_completed_forks();
      int N_active_forks = count_active_forks();
      
      std::cout << "\nprinting tree with "
	   << N_leaves
	   << " leaves and "
	   << N_forks
	   << " forks"
	   << std::endl;
      
      std::cout << "fork breakdown: "
	   << N_completed_forks
	   << " completed, "
	   << N_active_forks
	   << " active"
	   << std::endl;
		     
      std::cout << "leaves: "
	   << std::endl;
      
      for (std::string s: get_leaves())
	{
	  std::cout << s << " ";
	}
      std::cout << std::endl;

      if (N_completed_forks > 0)
	{
	  std::cout << "completed forks: "
	       << std::endl;
      
	  for (std::string s: get_completed_forks())
	    {
	      std::cout << s << " ";
	    }
	  std::cout << std::endl;
	}

      if (N_active_forks > 0)
	{
	  std::cout << "active forks: "
	       << std::endl;
      
	  for (std::string s: get_active_forks())
	    {
	      std::cout << s << " ";
	    }
	  std::cout << std::endl;
	}

      std::cout << "total_size = "
	   << total_size()
	   << std::endl;
      
      print_branch(root);
    }
  else
    {
      std::cout << "\n\ntree does not exist" << std::endl;
    }
  std::cout << "\n\n" << std::endl;
}

std::vector<std::string> btree::get_separable_forks()
{

  std::vector<std::string> separable_forks;

  // get the completed forks
  std::vector<std::string> comp_forks = get_completed_forks();

  std::cout << "completed forks" << std::endl;
  for (std::string fork : comp_forks)
    {
      std::cout << fork << std::endl;
    }


  // reduce the set of completed forks
  bool reduction_possible;
  std::string l_child, r_child;
  
  for (size_t i=0; i<comp_forks.size(); i++)
    {

      l_child = comp_forks[i] + "l";
      r_child = comp_forks[i] + "r";

      reduction_possible = false;
      for (size_t j=i+1; j<comp_forks.size(); j++)
	{
	  if (l_child == comp_forks[j])
	    {
	      reduction_possible = true;
	      break;
	    }
	}

      if (reduction_possible == false)
	{
	  separable_forks.push_back(l_child);
	}

      reduction_possible = false;
      for (size_t j=i+1; j<comp_forks.size(); j++)
	{
	  if (r_child == comp_forks[j])
	    {
	      reduction_possible = true;
	      break;
	    }
	}

      if (reduction_possible == false)
	{
	  separable_forks.push_back(r_child);
	}
      
    }

  if (separable_forks.size() == 0) separable_forks.push_back("m");

  return separable_forks;
  
}


int btree::get_generation(std::string loc)
{
  return get_branch(loc)->gen;
}


void btree::prepare_Nleaf_state(size_t N_leaves)
{
  int temp_size = single_size();
  size_t branch_idx;
  int min_gen, test_gen;

  destroy_tree();

  initialize_tree(temp_size);

  std::vector<std::string> current_leaves = get_leaves();

  while (current_leaves.size() < N_leaves)
    {

      for (size_t i=0; i<current_leaves.size(); i++)
	{
	  if (i == 0)
	    {
	      branch_idx = i;
	      min_gen = get_generation(current_leaves[i]);
	    }
	  else
	    {
	      test_gen = get_generation(current_leaves[i]);
	      if (test_gen < min_gen)
		{
		  branch_idx = i;
		  min_gen = test_gen;
		}
	    }
	}

      
      grow_at_branch_asym(current_leaves[branch_idx],
			  temp_size,
			  temp_size);
	
      current_leaves = get_leaves();
    }
  
}


void btree::foo()
{
  std::cout << "btree\n";
}
