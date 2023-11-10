#include <angle_array.hpp>

// constructor
angle_array::angle_array()
{
  N = -1;
  angles = nullptr;
}

// destructor
angle_array::~angle_array()
{
  destroy_angles();
}

// set the number of elements
void angle_array::set_N(int N)
{
  this->N = N;
  initialize_angles();
}

// get the number of elements
int angle_array::get_N()
{
  return N;
}

// initialize the array
void angle_array::initialize_angles()
{
  destroy_angles();

  if (N > 0)
    {

      angles = new angle[N];

      for (int i=0; i<N; i++)
	{
	  // angle id
	  angles[i].id = i + 1;
	  // type
	  angles[i].type = -1;
	  // angle atom ids
	  angles[i].i = -1;
	  angles[i].j = -1;
	  angles[i].k = -1;
	}
      
    }
  
}

// destroy the array
void angle_array::destroy_angles()
{
  if (angles != nullptr)
    {
      delete[] angles;
      angles = nullptr;
    }
}


// set angle element
void angle_array::set_angle(int i, angle a)
{
  angles[i] = a;
  angles[i].id = i + 1;
}


// get angle element
angle angle_array::get_angle(int i)
{
  return angles[i];
}


// write to stream
void angle_array::write(std::fstream &data_file)
{

  if (N > 0)
    {
      data_file << "\nAngles # angle-ID angle-type i j k\n" << std::endl;

      for (int i=0; i<N; i++)
	{
	  data_file << angles[i].id << "\t"
		    << angles[i].type << "\t"
		    << angles[i].i << "\t"
		    << angles[i].j << "\t"
		    << angles[i].k << std::endl;
	}
    }
  
}
