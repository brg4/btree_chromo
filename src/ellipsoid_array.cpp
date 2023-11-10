#include <ellipsoid_array.hpp>

// constructor
ellipsoid_array::ellipsoid_array()
{
  N = -1;
  ellipsoids = nullptr;
}

// destructor
ellipsoid_array::~ellipsoid_array()
{
  destroy_ellipsoids();
}

// set the number of elements
void ellipsoid_array::set_N(int N)
{
  if (N != this->N)
    {
      this->N = N;
      initialize_ellipsoids();
    }
}

// get the number of elements
int ellipsoid_array::get_N()
{
  return N;
}

// initialize the array
void ellipsoid_array::initialize_ellipsoids()
{
  destroy_ellipsoids();

  if (N > 0)
    {

      ellipsoids = new ellipsoid[N];

      for (int i=0; i<N; i++)
	{
	  // atom id
	  ellipsoids[i].id = i + 1;
	  // shape
	  ellipsoids[i].s.x = -1.0;
	  ellipsoids[i].s.y = -1.0;
	  ellipsoids[i].s.z = -1.0;
	  // orientation
	  ellipsoids[i].q.w = -1.0;
	  ellipsoids[i].q.v.x = -1.0;
	  ellipsoids[i].q.v.y = -1.0;
	  ellipsoids[i].q.v.z = -1.0;
	}
      
    }
  
}

// destroy the array
void ellipsoid_array::destroy_ellipsoids()
{
  if (ellipsoids != nullptr)
    {
      delete[] ellipsoids;
      ellipsoids = nullptr;
    }
}


// setter for elements
void ellipsoid_array::set_ellipsoid(int i, ellipsoid e)
{
  ellipsoids[i] = e;
}


// getter for elements
ellipsoid ellipsoid_array::get_ellipsoid(int i)
{
  return ellipsoids[i];
}


// set the minimum atom id and range accordingly [min,min+N)
void ellipsoid_array::set_min_id(int min_id)
{
  for (int i=0; i<N; i++)
    {
      ellipsoids[i].id = min_id + i;
    }
}


// set the shape of all ellipsoids
void ellipsoid_array::set_shape_all(vec shape)
{
  for (int i=0; i<N; i++)
    {
      ellipsoids[i].s = shape;
    }
}


// normalize all quaternions
void ellipsoid_array::normalize_quats()
{
  for (int i=0; i<N; i++)
    {
      ellipsoids[i].q = vqm.q_norm(ellipsoids[i].q);
    }
}


// set individual quaternion
void ellipsoid_array::set_quat(int i, quat q)
{
  ellipsoids[i].q.w = q.w;
  ellipsoids[i].q.v.x = q.v.x;
  ellipsoids[i].q.v.y = q.v.y;
  ellipsoids[i].q.v.z = q.v.z;
}


// set element-wise quaternions
void ellipsoid_array::set_quats(std::vector<quat> qs)
{
  for (int i=0; i<N; i++)
    {
      set_quat(i,qs[i]);
    }
}


// set the particle coordinates from a 1D array of doubles
void ellipsoid_array::set_quats_arr(double *&q, std::string order)
{

  quat p = vqm.q_null();
  
  if (order == "col")
    {
      for (int i=0; i<N; i++)
	{
	  p.w = q[i];
	  p.v.x = q[N+i];
	  p.v.y = q[2*N+i];
	  p.v.z = q[3*N+i];
	  set_quat(i,p);
	}
    }
  else if (order == "row")
    {
      for (int i=0; i<N; i++)
	{
	  p.w = q[4*i];
	  p.v.z = q[4*i+1];
	  p.v.y = q[4*i+2];
	  p.v.z = q[4*i+3];
	  set_quat(i,p);
	}
    }
}


// read coordinates from binary file
int ellipsoid_array::read_bin_quats(std::string data_filename, std::string order, bool force_resize)
{
  std::fstream data_file;

  data_file.open(data_filename, std::ios::in | std::ios::binary | std::ios::ate);

  // int data_size = 4;
    
  if (!data_file.is_open())
    {
      std::cout << "ERROR: file not opened in read_bin_quats" << std::endl;
      return 1;
    }
  else
    {

      // read in the binary file
      std::streampos size = data_file.tellg();
      char *memblock;
      double *q;

      int data_size = sizeof(double);
      char vals[sizeof(double)];
      int N_data = size/data_size;
      std::cout << "N_data = " << N_data << std::endl;
      int N_bin_quats = N_data/4;
      std::cout << "N_bin_quats = " << N_bin_quats << std::endl;

      // test for resizing
      if (force_resize)
	{
	  set_N(N_bin_quats);
	}
      else
	{
	  if (N_bin_quats != get_N())
	    {
	      std::cout << "ERROR: incorrect ellipsoid_array size" << std::endl;
	      return 1;
	    }
	}

      memblock = new char[size];
      
      data_file.seekg(0, std::ios::beg);

      data_file.read(memblock,size);

      data_file.close();

      q = new double[N_data];

      for (int i=0; i<N_data; i++)
	{
	  for (int j=0; j<data_size; j++)
	    {
	      vals[j] = memblock[i*data_size+j];
	    }
	  memcpy(&q[i], vals, 8);
	}

      delete[] memblock;
      
      set_quats_arr(q,order);

      delete[] q;
      return 0;
      
    }
    
}


// write to stream
void ellipsoid_array::write(std::fstream &data_file)
{

  if (N > 0)
    {
      data_file << "\nEllipsoids # atom-ID shapex shapey shapez quatw quati quatj quatk\n" << std::endl;

      for (int i=0; i<N; i++)
	{        
	  data_file << ellipsoids[i].id << "\t"
		    << ellipsoids[i].s.x << "\t"
		    << ellipsoids[i].s.y << "\t"
		    << ellipsoids[i].s.z << "\t"
		    << ellipsoids[i].q.w << "\t"
		    << ellipsoids[i].q.v.x << "\t"
		    << ellipsoids[i].q.v.y << "\t"
		    << ellipsoids[i].q.v.z << std::endl;
	}
    }
  
}


// write quaternions to a binary file
int ellipsoid_array::write_bin(std::string data_filename, std::string order)
{
  std::fstream data_file;

  data_file.open(data_filename, std::ios::out | std::ios::binary);

  // int data_size = 4;
    
  if (!data_file.is_open())
    {
      std::cout << "ERROR: file not opened in write_bin" << std::endl;
      return 1;
    }
  else if (N > 0)
    {
      double *q = new double[4*N];

      if (order == "col")
	{
	  for (int i=0; i<N; i++)
	    {
	      q[i] = ellipsoids[i].q.w;
	      q[N+i] = ellipsoids[i].q.v.x;
	      q[2*N+i] = ellipsoids[i].q.v.y;
	      q[3*N+i] = ellipsoids[i].q.v.z;
	    }
	}
      else if (order == "row")
	{
	  for (int i=0; i<N; i++)
	    {
	      q[4*i] = ellipsoids[i].q.w;
	      q[4*i+1] = ellipsoids[i].q.v.x;
	      q[4*i+2] = ellipsoids[i].q.v.y;
	      q[4*i+3] = ellipsoids[i].q.v.z;
	    }
	}

      data_file.write(reinterpret_cast<const char*>(q),4*N*sizeof(double));
      
      delete[] q;
      data_file.close();
      return 0;
    }
  else
    {
      std::cout << "no quats to write" << std::endl;
      data_file.close();
      return 0;
    }

  return 0;


}
