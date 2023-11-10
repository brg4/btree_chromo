#include <loop_topology.hpp>

// constructor
loop_topology::loop_topology()
{
  rand_eng.seed(0);
}


// destructor
loop_topology::~loop_topology()
{
}


void loop_topology::prng_seed(int s)
{
  rand_eng.seed(s);
}


void loop_topology::set_step_dist(std::string family, double l, int k_max)
{
  if (family == "poisson")
    {
      step_dist.poisson_dist(l,k_max+1);
    }
  else if (family == "uniform")
    {
      step_dist.uniform_dist(k_max+1);
    }
}


void loop_topology::set_coords(std::vector<vec> coords)
{
  this->coords = coords;
}


// prepare the vector of binding regions
void loop_topology::prepare_binding_regions(std::vector<std::string> leaves, std::vector<theta_topo> leaf_topos, int *&t)
{

  regions.clear();

  std::vector<int> partitions;

  std::string leaf; // leaf that binding region belongs to
  int ll, ul, size; // lower limit and upper limit of indices
  bool completed;
  bool ter_crossing;
  int mid_ll = -1;
  int mid_ul = -1;

  theta_topo leaf_topo;

  for (size_t i_leaf=0; i_leaf<leaves.size(); i_leaf++)
    {

      leaf = leaves[i_leaf];
      leaf_topo = leaf_topos[i_leaf];

      // test if the leaf is a completed chromosome
      if ((leaf_topo.start == leaf_topo.end_link) &&
	  (leaf_topo.end == leaf_topo.start_link))
	{

	  // determine the partitions along the leaf
	  partitions.clear();

	  for (int i=(leaf_topo.start+1); i<leaf_topo.end; i++)
	    {
	      if (t[i] == 6)
		{
		  partitions.push_back(i);
		}
	    }

	  // prepare a parititioning if forks are present
	  if (partitions.size() > 0)
	    {
	    
	      // create regions from the partitions
	      for (size_t i=0; i<(partitions.size()-1); i++)
		{
		  if ((partitions[i+1] - partitions[i]) > 3)
		    {
		      ll = partitions[i] + 1;
		      ul = partitions[i+1] - 1;
		      size = ul - ll + 1;
		      completed = false;
		      ter_crossing = false;
		      mid_ll = -1;
		      mid_ul = -1;
		      binding_region region(leaf,ll,ul,size,
					    completed,
					    ter_crossing,mid_ll,mid_ul);
		      regions.push_back(region);
		    }
		}

	      // test for crossing back over Ter
	      int completed_gap = 0;

	      completed_gap += (leaf_topo.end - partitions[partitions.size()-1]);
	      completed_gap += (partitions[0] - leaf_topo.start);

	      if (completed_gap > 1)
		{
		  ter_crossing = true;
		  if (partitions[partitions.size()-1] < leaf_topo.end)
		    {
		      ul = partitions[partitions.size()-1] + 1;
		      mid_ul = leaf_topo.end;
		    }
		  else
		    {
		      ter_crossing = false;
		      ll = leaf_topo.start;
		      ul = partitions[0] - 1;
		    }
		  if (partitions[0] > leaf_topo.start)
		    {
		      ll = partitions[0] - 1;
		      mid_ll = leaf_topo.start;
		    }
		  else
		    {
		      ter_crossing = false;
		      ul = leaf_topo.end;
		      ll = partitions[partitions.size()-1] + 1;
		    }
		  completed = false;
		  size = completed_gap;
		  binding_region region(leaf,ll,ul,size,
					completed,
					ter_crossing,mid_ll,mid_ul);
		  regions.push_back(region);
		}
	      
	    }
	  else // no forks are present
	    {

	      ll = leaf_topo.start;
	      ul = leaf_topo.end;
	      size = ul - ll + 1;
	      completed = true;
	      ter_crossing = false;
	      mid_ll = -1;
	      mid_ul = -1;
	      binding_region region(leaf,ll,ul,size,
				    completed,
				    ter_crossing,mid_ll,mid_ul);
	      regions.push_back(region);
	      
	    }
	  
	}
      else
	{

	  // determine the partitions along the leaf
	  partitions.clear();

	  partitions.push_back(leaf_topo.start);
	  
	  for (int i=(leaf_topo.start+1); i<leaf_topo.end; i++)
	    {
	      if (t[i] == 6)
		{
		  partitions.push_back(i);
		}
	    }

	  partitions.push_back(leaf_topo.end);

	  // create regions from the partitions
	  for (size_t i=0; i<(partitions.size()-1); i++)
	    {
	      if ((partitions[i+1] - partitions[i]) > 3)
		{
		  ll = partitions[i] + 1;
		  ul = partitions[i+1] - 1;
		  size = ul - ll + 1;
		  completed = false;
		  ter_crossing = false;
		  mid_ll = -1;
		  mid_ul = -1;
		  binding_region region(leaf,ll,ul,size,
					completed,
					ter_crossing,mid_ll,mid_ul);
		  regions.push_back(region);
		}
	    }
	  
	} // end conditinal for complete leaf test
      
    } // end loop over leaves

  //std::cout << "at end of prepare_binding_regions" << std::endl;
  for (size_t i_reg=0; i_reg<regions.size(); i_reg++)
    {
      //std::cout << "i_reg = " << i_reg << std::endl;
      regions[i_reg].prepare_idx();
      //regions[i_reg].print_region_map();
    }
  
}


// initialize the loops based on the set of binding regions
void loop_topology::initialize_loops(int N_loops, int min_dist)
{
  // clear the loops and binding regions
  loops.clear();
  
  size_t N_regions = regions.size();

  // std::cout << "at start of initialize_loops" << std::endl;
  // for (size_t i_reg=0; i_reg<regions.size(); i_reg++)
  //   {
  //     std::cout << "i_reg = " << i_reg << std::endl;
  //     regions[i_reg].print_region_map();
  //   }

  // only select binding regions that can contain a full loop (hinge and anchor) upon initialization

  int total_binding_region_size = 0;

  for (size_t i_region=0; i_region<N_regions; i_region++)
    {
      if (regions[i_region].get_size() > 2*min_dist)
	{
	  total_binding_region_size += regions[i_region].get_size();
	}
    }

  std::uniform_int_distribution<int> unif_dist(1,total_binding_region_size);
  std::vector<int> a_dist;

  // loop over the number of loops
  for (int i_loop=0; i_loop<N_loops; i_loop++)
    {
      a_dist.push_back(unif_dist(rand_eng));
    }

  // select an anchor region
  for (int i_loop=0; i_loop<N_loops; i_loop++)
    {
      int accumulator = 0;
      for (size_t i_region=0; i_region<N_regions; i_region++)
	{
	  if (regions[i_region].get_size() > 2*min_dist)
	    {
	      accumulator += regions[i_region].get_size();
	    } 
	  if (a_dist[i_loop] <= accumulator)
	    {
	      loop l(i_region,i_region);
	      loops.push_back(l);
	      break;
	    }
	}
    }


  int a_reg, h_reg;

  std::uniform_real_distribution<double> dir_dist(0.0, 1.0);
  double r_d;
  int d, r_mono, a_mono, h_mono;
  
  // select an anchor for each loop
  for (int i_loop=0; i_loop<N_loops; i_loop++)
    {
      // get the anchor region
      a_reg = loops[i_loop].get_a_region();

      // select a random monomer for the anchor
      std::uniform_int_distribution<int> unif_dist(0,regions[a_reg].get_size()-1);
      r_mono = unif_dist(rand_eng);
      a_mono = regions[a_reg].get_mono_pos(r_mono);
      loops[i_loop].set_a(a_mono);
      
      // bind the anchor
      loops[i_loop].set_a_bound(true);

      // get the hinge region
      h_reg = loops[i_loop].get_h_region();

      // based on the position of the anchor and a minimum distance between the anchor and hinge, select a direction for the hinge to travel
      r_d = dir_dist(rand_eng);
      d = regions[h_reg].select_direction(a_mono,
					  min_dist,
					  r_d);
      loops[i_loop].set_d(d);

      // select a hinge based on the anchor position and direction
      h_mono = regions[h_reg].get_relative_monomer_pos(a_mono,
						       min_dist,
						       d);
      loops[i_loop].set_h(h_mono);

      // std::cout << "a=" << a_mono << ", h=" << h_mono << ", d=" << d << std::endl;
      // regions[a_reg].print_region_map();
      
      // bind the hinge
      loops[i_loop].set_h_bound(true);
    }

}


// prepare the vector of binding regions
void loop_topology::update_loops(int ext_max, int min_dist, double p_unbinding, double r_g)
{

  std::uniform_real_distribution<double> unif_dist(0.0, 1.0);
  std::vector<int> intra_updates;
  std::vector<std::vector<int>> inter_updates;
  int a_mono, h_mono, d;
  size_t a_reg, h_reg;
  vec a_coord;
  bool h_bound;
  double r_intra, r_d;
  size_t r_inter, accumulator;
  int h_intra, h_inter;
  
  // loop over the loops
  for (size_t i_loop=0; i_loop<loops.size(); i_loop++)
    {

      // get the anchor, anchor region, and anchor coordinate
      a_mono = loops[i_loop].get_a();
      a_reg = static_cast<size_t>(loops[i_loop].get_a_region());
      a_coord = coords[a_mono];

      // get the hinge, hinge region, and hinge binding state
      h_mono = loops[i_loop].get_h();
      h_reg = static_cast<size_t>(loops[i_loop].get_h_region());
      h_bound = loops[i_loop].get_h_bound();

      // hinge is already bound
      if (h_bound == true)
	{

	  // get the direction
	  d = loops[i_loop].get_d();

	  // update the proximities in the hinge region
	  // std::cout << "h_reg = " << h_reg << std::endl;
	  // regions[h_reg].print_region_map();
	  regions[h_reg].update_proximities(r_g,a_coord,coords);

	  if (a_reg == h_reg) regions[h_reg].filter_proximities_near_a(min_dist,a_mono);
	  
	  // get the intra updates
	  intra_updates = regions[h_reg].get_and_filter_intra_candidates(ext_max,h_mono,d);

	  // std::cout << "intra_updates, " << h_mono << ", " << intra_updates.size() << " candidates" << std::endl;
	  // for (size_t i_intra=0; i_intra<intra_updates.size(); i_intra++)
	  //   {
	  //     std::cout << intra_updates[i_intra] << std::endl;
	  //   }

	  if (unif_dist(rand_eng) > p_unbinding)
	    {
	  
	      // select an updated hinge for 1D motion along strand
	      r_intra = unif_dist(rand_eng);
	      std::cout << "r_intra = " << r_intra << std::endl;
	      h_intra = intra_updates[step_dist.get_k(r_intra,static_cast<int>(intra_updates.size()))];
	      std::cout << "h_intra_final = " << h_intra << std::endl;

	      loops[i_loop].set_h(h_intra);
	  
	    }
	  else
	    {
	  
	      // select an updated hinge for 3D motion between strands

	      int N_inter_total = 0;

	      // determine the candidate inter-strand hinge updates for all the regions
	      for (size_t i_reg=0; i_reg<regions.size(); i_reg++)
		{
		  // update proximities not in the hinge region
		  if (i_reg != h_reg)
		    {
		      regions[i_reg].update_proximities(r_g,a_coord,coords);
		    }

		  // filter proximities violating minimum distance within anchor region
		  if (i_reg == a_reg)
		    {
		      regions[i_reg].filter_proximities_near_a(min_dist,a_mono);
		    }

		  // add the region's candidates to the total set of possible inter updates
		  inter_updates.push_back(regions[i_reg].get_inter_candidates());

		  N_inter_total += inter_updates[i_reg].size();
	      
		} // end loop to determine inter candidates

	      if (N_inter_total > 0)
		{
		  std::uniform_int_distribution<size_t> unif_dist(0,N_inter_total);

		  accumulator = 0;
		  r_inter = unif_dist(rand_eng);
		  std::cout << "r_inter = " << r_inter  << ", N_inter_total = " << N_inter_total << std::endl;

		  // determine the candidate inter-strand hinge updates for all the regions
		  for (size_t i_reg=0; i_reg<regions.size(); i_reg++)
		    {

		      if ((accumulator + inter_updates[i_reg].size()) > r_inter)
			{
			  h_inter = static_cast<int>(r_inter-accumulator);
			  // set the hinge
			  std::cout << "accumulator = " << accumulator  << ", i_reg = " << i_reg << std::endl;
			  std::cout << "size = " << inter_updates[i_reg].size()  << ", h_inter = " << h_inter << std::endl;
			  loops[i_loop].set_h(inter_updates[i_reg][h_inter]);
			  // set the hinge region
			  loops[i_loop].set_h_region(i_reg);
			  // bind the hinge
			  loops[i_loop].set_h_bound(true);
			  // set the direction
			  d = 1;
			  r_d = unif_dist(rand_eng);
			  if (r_d < 0.5) d = -1;
			  loops[i_loop].set_d(d);
			}
		      else
			{
			  accumulator += inter_updates[i_reg].size();
			}
	      
		    }
		  
		  
		}
	      else
		{

		  // make the hinge unbound
		  loops[i_loop].set_h_bound(false);
		  
		}

	    }
	  
	}
      else // hinge is unbound
	{

	  // select an updated hinge for 3D motion between strands
	  int N_inter_total = 0;

	  // determine the candidate inter-strand hinge updates for all the regions
	  for (size_t i_reg=0; i_reg<regions.size(); i_reg++)
	    {
	      // update proximities for all regions
	      regions[i_reg].update_proximities(r_g,a_coord,coords);

	      // filter proximities violating minimum distance within anchor region
	      if (i_reg == a_reg)
		{
		  regions[i_reg].filter_proximities_near_a(min_dist,a_mono);
		}

	      // add the region's candidates to the total set of possible inter updates
	      inter_updates.push_back(regions[i_reg].get_inter_candidates());

	      N_inter_total += inter_updates[i_reg].size();
	      
	    } // end loop over regions

	  if (N_inter_total > 0)
	    {
	      std::uniform_int_distribution<int> unif_dist(0,N_inter_total);

	      accumulator = 0;
	      r_inter = unif_dist(rand_eng);

	      // determine the candidate inter-strand hinge updates for all the regions
	      for (size_t i_reg=0; i_reg<regions.size(); i_reg++)
		{

		  if ((accumulator + inter_updates[i_reg].size()) >= r_inter)
		    {
		      h_inter = static_cast<int>(r_inter-accumulator);
		      // set the hinge
		      loops[i_loop].set_h(inter_updates[i_reg][h_inter]);
		      // set the hinge region
		      loops[i_loop].set_h_region(i_reg);
		      // bind the hinge
		      loops[i_loop].set_h_bound(true);
		      // set the direction
		      d = 1;
		      r_d = unif_dist(rand_eng);
		      if (r_d < 0.5) d = -1;
		      loops[i_loop].set_d(d);
		    }
		  else
		    {
		      accumulator += inter_updates[i_reg].size();
		    }
	      
		}
	      
	    } // end conditional for nonzero inter candidates
	  
	} // end conditional for bound/unbound hinges
      
    } // end loop over loops
  
}


// get the loops
std::vector<loop> loop_topology::get_loops()
{
  return loops;
}

// get the binding regions
std::vector<binding_region> loop_topology::get_regions()
{
  return regions;
}

