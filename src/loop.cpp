#include <loop.hpp>

// constructor
loop::loop(int a_region, int h_region)
{
  set_a_region(a_region);
  set_h_region(h_region);
  set_a_bound(false);
  set_h_bound(false);
  set_a(0);
  set_h(0);
  set_d(0);
}


// destructor
loop::~loop()
{
}


// set a binding state
void loop::set_a_bound(bool a_bound)
{
  this->a_bound = a_bound;
}


// set h binding state
void loop::set_h_bound(bool h_bound)
{
  this->h_bound = h_bound;
}


// set a
void loop::set_a(int a)
{
  this->a = a;
}


// set h
void loop::set_h(int h)
{
  this->h = h;
}


// set d
void loop::set_d(int d)
{
  this->d = d;
}


// set a_region
void loop::set_a_region(int a_region)
{
  this->a_region = a_region;
}


// set h_region
void loop::set_h_region(int h_region)
{
  this->h_region = h_region;
}


// get a binding state
bool loop::get_a_bound()
{
  return a_bound;
}


// get a binding state
bool loop::get_h_bound()
{
  return h_bound;
}


// get a
int loop::get_a()
{
  return a;
}


// get h
int loop::get_h()
{
  return h;
}


// get d
int loop::get_d()
{
  return d;
}


// get a_region
int loop::get_a_region()
{
  return a_region;
}


// get h_region
int loop::get_h_region()
{
  return h_region;
}
