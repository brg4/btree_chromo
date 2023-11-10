#include <boundary_surface.hpp>

// constructor
boundary_surface::boundary_surface()
{
}

// destructor
boundary_surface::~boundary_surface()
{
}


// create a new tri_face
tri_face boundary_surface::new_tri_face(int v0, int v1, int v2)
{
  tri_face f;
  f.verts[0] = v0;
  f.verts[1] = v1;
  f.verts[2] = v2;
  return f;
}


// project the coords onto the surface of a sphere
void boundary_surface::project_to_sphere()
{
  for (size_t i=0; i<coords.size(); i++)
    {
      coords[i] = vqm.v_norm(coords[i]);
    }
}


// project the coords onto the surface of a sphere
void boundary_surface::scale_coords(double s)
{
  for (size_t i=0; i<coords.size(); i++)
    {
      coords[i] = vqm.v_ax(s,coords[i]);
    }
}


// initialize a unit icosahedron
void boundary_surface::unit_icosahedron()
{
  const double r[3] = {0.0,1.0,(1.0+sqrt(5.0))/2.0};
  const int x_permut[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
  const int s_permut[4][3] = {{0,1,1},{0,-1,1},{0,-1,-1},{0,1,-1}};

  tri_surf.clear();
  coords.clear();

  vec temp_r;

  for (int i=0; i<3; i++)
    {
      for (int j=0; j<4; j++)
	{

	  temp_r = vqm.v_new(s_permut[j][x_permut[i][0]]*r[x_permut[i][0]],
			     s_permut[j][x_permut[i][1]]*r[x_permut[i][1]],
			     s_permut[j][x_permut[i][2]]*r[x_permut[i][2]]);

	  coords.push_back(temp_r);

	}
    }


  // upper pole
  tri_surf.push_back(new_tri_face(0,1,8));
  tri_surf.push_back(new_tri_face(0,8,4));
  tri_surf.push_back(new_tri_face(0,4,5));
  tri_surf.push_back(new_tri_face(0,5,11));
  tri_surf.push_back(new_tri_face(0,11,1));

  // lower pole
  tri_surf.push_back(new_tri_face(2,3,9));
  tri_surf.push_back(new_tri_face(2,9,7));
  tri_surf.push_back(new_tri_face(2,7,6));
  tri_surf.push_back(new_tri_face(2,6,10));
  tri_surf.push_back(new_tri_face(2,10,3));

  // upper equator
  tri_surf.push_back(new_tri_face(1,7,8));
  tri_surf.push_back(new_tri_face(8,9,4));
  tri_surf.push_back(new_tri_face(4,3,5));
  tri_surf.push_back(new_tri_face(5,10,11));
  tri_surf.push_back(new_tri_face(11,6,1));

  // lower equator
  tri_surf.push_back(new_tri_face(3,4,9));
  tri_surf.push_back(new_tri_face(9,8,7));
  tri_surf.push_back(new_tri_face(7,1,6));
  tri_surf.push_back(new_tri_face(6,11,10));
  tri_surf.push_back(new_tri_face(10,5,3));
  

  project_to_sphere();
  
}


// initialize a unit tetrahedron
void boundary_surface::unit_tetrahedron()
{
  tri_surf.clear();
  coords.clear();


  coords.push_back(vqm.v_new(sqrt(8.0/9.0),0.0,-1.0/3.0));
  coords.push_back(vqm.v_new(-sqrt(2.0/9.0),sqrt(2.0/3.0),-1.0/3.0));
  coords.push_back(vqm.v_new(-sqrt(2.0/9.0),-sqrt(2.0/3.0),-1.0/3.0));
  coords.push_back(vqm.v_new(0.0,0.0,1.0));

  tri_surf.push_back(new_tri_face(0,1,2));
  tri_surf.push_back(new_tri_face(1,2,3));
  tri_surf.push_back(new_tri_face(2,3,0));
  tri_surf.push_back(new_tri_face(0,1,3));

}


// interpolate triangles within a single face
void boundary_surface::interpolate_face(tri_face t_f, std::vector<edge_map> &unique_edges)
{
  
  const int perm_edges[3][2] = {{0,1},{1,2},{2,0}};
  int new_verts[3];
  std::array<int,2> temp_edge;

  for (int i=0; i<3; i++)
    {
      temp_edge[0] = t_f.verts[perm_edges[i][0]];
      temp_edge[1] = t_f.verts[perm_edges[i][1]];

      new_verts[i] = vert_from_edge(temp_edge,unique_edges);

      if (new_verts[i] == -1) std::cout << temp_edge[0] << "," << temp_edge[1] << " fail" << std::endl;
    }

  //   0
  //  3 5
  // 1 4 2

  // upper
  tri_surf.push_back(new_tri_face(t_f.verts[0],new_verts[0],new_verts[2]));
  // left
  tri_surf.push_back(new_tri_face(t_f.verts[1],new_verts[1],new_verts[0]));
  // right
  tri_surf.push_back(new_tri_face(t_f.verts[2],new_verts[2],new_verts[1]));
  // center triangle
  tri_surf.push_back(new_tri_face(new_verts[0],new_verts[1],new_verts[2]));

}


// test edge equivalence
bool boundary_surface::edge_equiv(std::array<int,2> &e0, std::array<int,2> &e1)
{
  if (((e0[0] == e1[0]) && (e0[1] == e1[1])) ||
      ((e0[1] == e1[0]) && (e0[0] == e1[1])))
    {
      return true;
    }
  return false;
}


// find vertex from edge mapping
int boundary_surface::vert_from_edge(std::array<int,2> &e, std::vector<edge_map> &edge_mapping)
{
  for (edge_map e_m : edge_mapping)
    {
      if (edge_equiv(e,e_m.edge)) return e_m.vert;
    }
  return -1;
}


// interpolate along the current surfaces
void boundary_surface::interpolate_surface()
{

  const int perm_edges[3][2] = {{0,1},{1,2},{2,0}};
  std::vector<tri_face> old_tri_surf = tri_surf;
  tri_surf.clear();

  // determine the set of unique edges between vertices
  std::vector<std::array<int,2>> unique_edges;
  bool new_edge;
  std::array<int,2> temp_edge;

  unique_edges.clear();

  // loop over old faces
  for (tri_face t_f : old_tri_surf)
    {

      // std::cout << t_f.verts[0] << "," << t_f.verts[1] << "," << t_f.verts[2] << std::endl;
      // loop over edges in face of old surface
      for (int j=0; j<3; j++)
	{
	  
	  temp_edge[0] = t_f.verts[perm_edges[j][0]];
	  temp_edge[1] = t_f.verts[perm_edges[j][1]];

	  // loop over set of unique edges
	  new_edge = true;
	  for (std::array<int,2> edge : unique_edges)
	    {
	      if (edge_equiv(temp_edge,edge) == true)
		{
		  new_edge = false;
		  break;
		}
	    }

	  // add edge if unique
	  if (new_edge == true)
	    {
	      unique_edges.push_back(temp_edge);
	      // std::cout << temp_edge[0] << " " << temp_edge[1] << std::endl;
	    }
	  
	}
    }
  edge_map temp_e_m;
  std::vector<edge_map> edge_mapping;

  edge_mapping.clear();
  
  for (std::array<int,2> edge : unique_edges)
    {
      // std::cout << edge[0] << "," << edge[1] << std::endl;
      coords.push_back((vqm.v_linterp(0.5,
				      coords[edge[0]],
				      coords[edge[1]])));
      temp_e_m.vert = coords.size() - 1;
      temp_e_m.edge = edge;
      edge_mapping.push_back(temp_e_m);
    }

  // std::cout << "number old_tri_surfs = " << old_tri_surf.size() << std::endl;
  
  for (tri_face t_f : old_tri_surf)
    {
      // std::cout << t_f.verts[0] << "," << t_f.verts[1] << "," << t_f.verts[2] << std::endl;
      interpolate_face(t_f,edge_mapping);
      // std::cout << tri_surf.size() << std::endl;
    }

  // std::cout << "number tri_surfs = " << tri_surf.size() << std::endl;
  
}


void boundary_surface::generate_sphere(double R, double r)
{

  size_t N_target = 4*pow((R+r)/r,2.0);
  
  unit_icosahedron();

  while (coords.size() < N_target)
    {
      interpolate_surface();
      project_to_sphere();
    }
  scale_coords(R+r);
  
}


// getter for coordinates
std::vector<vec> boundary_surface::get_coords()
{
  return coords;
}


// getter for number of vertices
int boundary_surface::get_N_verts()
{
  return static_cast<int>(coords.size());
}


// write the boundary coordinates to an xyz file
void boundary_surface::write_xyz(std::string data_filename)
{

  // begin writing data file
  
  std::fstream data_file;

  data_file.open(data_filename, std::ios::out);

  if (!data_file)
    {
      std::cout << "ERROR: file not opened in write_xyz" << std::endl;
    }
  else
    {

      // write system summary
      data_file << coords.size() << "\n" << std::endl;

      for (vec r : coords)
	{
	  data_file << "C\t"
		    << r.x << "\t"
		    << r.y << "\t"
		    << r.z << std::endl;
	}
      
    }

  data_file.close();
  
}
