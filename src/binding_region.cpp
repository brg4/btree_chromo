#include <binding_region.hpp>

// constructor
binding_region::binding_region(std::string leaf, int ll, int ul, int size,
		 bool completed,
		 bool ter_crossing, int mid_ll, int mid_ul)
{
  this->leaf = leaf;
  this->ll = ll;
  this->ul = ul;
  this->size = size;
  this->completed = completed;
  this->ter_crossing = ter_crossing;
  this->mid_ll = mid_ll;
  this->mid_ul = mid_ul;
  
}


// destructor
binding_region::~binding_region()
{
  if (mono_idx != nullptr)
    {
      delete[] mono_idx;
      mono_idx = nullptr;
    }
  if (proximities != nullptr)
    {
      delete[] proximities;
      proximities = nullptr;
    }
  reg_idx.clear();
}


// prepare the array and map for monomer and region indices
void binding_region::prepare_idx()
{
  if (mono_idx != nullptr)
    {
      delete[] mono_idx;
      mono_idx = nullptr;
    }
  
  mono_idx = new int[size];

  reg_idx.clear();

  if (ter_crossing == true)
    {
      int i = 0;
      for (int j=ul; j<(mid_ul+1); j++)
	{
	  mono_idx[i] = j;
	  reg_idx[j] = i;
	  i += 1;
	}
      for (int j=mid_ll; j<(ll+1); j++)
	{
	  mono_idx[i] = j;
	  reg_idx[j] = i;
	  i += 1;
	}
    }
  else
    {
      for (int i=0; i<size; i++)
	{
	  mono_idx[i] = ll + i;
	  reg_idx[ll+i] = i;
	}
    }

  if (proximities != nullptr)
    {
      delete[] proximities;
      proximities = nullptr;
    }

  proximities = new int[size];

  reset_proximities();
}


void binding_region::print_region_map()
{
  std::cout << "size = " << size << std::endl;
  int j;
  for (int i=0; i<size; i++)
    {
      j = get_mono_pos(i);
      std::cout << i << "," << j << "," << get_reg_pos(j) << std::endl;
    }
}


// getter for size
int binding_region::get_size()
{
  return size;
}


// getter for leaf
std::string binding_region::get_leaf()
{
  return leaf;
}


// getter for region position from monomer position
int binding_region::get_reg_pos(int mono_pos)
{
  return reg_idx[mono_pos];
}


// getter for monomer position from region position
int binding_region::get_mono_pos(int reg_pos)
{
  return mono_idx[reg_pos];
}


// based on the position of the anchor and a minimum distance between the anchor and hinge, select a direction for the hinge to travel
int binding_region::select_direction(int mono_pos, int min_dist, double r_d)
{
  // randomly sample a point within the entire region
  int d;
  
  if (completed == true)
    {
      d = 1;
      if (r_d < 0.5) d = -1;
    }
  else
    {
      int reg_pos = get_reg_pos(mono_pos);

      if (reg_pos < min_dist)
	{
	  d = 1;
	}
      else if (reg_pos >= (size - min_dist))
	{
	  d = -1;
	}
      else
	{
	  d = 1;
	  if (r_d < 0.5) d = -1;
	}
    }
  return d;
}


// determine a new relative pos based on the current position, distance, and direction
int binding_region::get_relative_monomer_pos(int mono_pos, int dist, int dir)
{

  if (dist == 0) return mono_pos;
  
  int reg_pos = get_reg_pos(mono_pos) + dir*dist;
  
  if (completed == true)
    {
      if (reg_pos < 0)
	{
	  reg_pos = size + reg_pos;
	}
      else if (reg_pos >= size)
	{
	  reg_pos = reg_pos - size;
	}
    }
  else
    {
      if (reg_pos < 0)
	{
	  reg_pos = 0;
	}
      else if (reg_pos >= size)
	{
	  reg_pos = size - 1;
	}
    }
  return get_mono_pos(reg_pos);
}


// get the distance between two monomer positions in a region
int binding_region::get_dist(int mono_pos_i, int mono_pos_j)
{
  int reg_pos_i = get_reg_pos(mono_pos_i);
  int reg_pos_j = get_reg_pos(mono_pos_j);
  int d;

  if (completed == true)
    {
      if (reg_pos_i > reg_pos_j)
	{
	  d = std::min(reg_pos_i-reg_pos_j,size-reg_pos_i+reg_pos_j);
	}
      else if (reg_pos_j > reg_pos_i)
	{
	  d = std::min(reg_pos_j-reg_pos_i,size-reg_pos_j+reg_pos_i);
	}
      else
	{
	  d = 0;
	}
    }
  else
    {
      d = abs(reg_pos_i - reg_pos_j);
    }
  return d;
}


// reset the proximities
void binding_region::reset_proximities()
{
  for (int i=0; i<size; i++)
    {
      proximities[i] = 0;
    }
}


// update the proximities given a grab radius and coordinate
void binding_region::update_proximities(double r_g, vec a_coord, std::vector<vec> &coords)
{
  double r_g_2 = pow(r_g,2.0);
  double test_dist_L2;
  
  for (int i=0; i<size; i++)
    {
      // std::cout << "i = " << i << std::endl;
      // std::cout << "m_i = " << get_mono_pos(i) << std::endl;
      test_dist_L2 = vqm.v_L2(vqm.v_xpy(a_coord,vqm.v_inv(coords[get_mono_pos(i)])));
      if (test_dist_L2 < r_g_2)
	{
	  proximities[i] = 1;
	}
      else
	{
	  proximities[i] = 0;
	}
    }
}


// filter the proximites too close to the anchor
void binding_region::filter_proximities_near_a(int min_dist, int a_mono_pos)
{
  int i_reg_pos;
  i_reg_pos = get_reg_pos(a_mono_pos);
  proximities[i_reg_pos] = 0;

  for (int i=1; i<(min_dist+1); i++)
    {
      // reverse direction
      i_reg_pos = get_reg_pos(get_relative_monomer_pos(a_mono_pos,i,-1));
      proximities[i_reg_pos] = 0;
      // forward direction
      i_reg_pos = get_reg_pos(get_relative_monomer_pos(a_mono_pos,i,1));
      proximities[i_reg_pos] = 0;
    }
}


// get the intra candidates and remove them from the set of proximities
std::vector<int> binding_region::get_and_filter_intra_candidates(int ext_max, int h_mono_pos, int dir)
{
  std::vector<int> intra_candidates;
  intra_candidates.push_back(h_mono_pos);
  proximities[get_reg_pos(h_mono_pos)] = 0;

  int i_reg_pos, i_mono_pos;
  int prev_mono_pos = h_mono_pos;


  for (int i=1; i<(ext_max+1); i++)
    {
      i_mono_pos = get_relative_monomer_pos(h_mono_pos,i,dir);
      i_reg_pos = get_reg_pos(i_mono_pos);

      if ((proximities[i_reg_pos] == 1) && (i_mono_pos != prev_mono_pos))
	{
	  intra_candidates.push_back(i_mono_pos);
	  prev_mono_pos = i_mono_pos;
	  proximities[i_reg_pos] = 0;
	}
      else
	{
	  break;
	}
    }
      
  return intra_candidates;
}

// get the inter candidates
std::vector<int> binding_region::get_inter_candidates()
{
  std::vector<int> inter_candidates;

  for (int i=0; i<size; i++)
    {
      if (proximities[i] == 1) inter_candidates.push_back(get_mono_pos(i));
    }
  
  return inter_candidates;
}
