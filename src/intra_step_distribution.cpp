#include <intra_step_distribution.hpp>

// constructor
intra_step_distribution::intra_step_distribution()
{
  N = 0;
  w = nullptr;
  cw = nullptr;
}


// destructor
intra_step_distribution::~intra_step_distribution()
{
  destroy_w_cw();
}


// size the distribution
void intra_step_distribution::set_N(int N)
{
  this->N = N;
  reset_w_cw();
}


// reset the weights
void intra_step_distribution::reset_w_cw()
{
  destroy_w_cw();
  w = new double[N];
  cw = new double[N];
  for (int i=0; i<N; i++)
    {
      w[i] = 0.0;
      cw[i] = 0.0;
    }
}


// destroy the weight array
void intra_step_distribution::destroy_w_cw()
{
  if (w != nullptr)
    {
      delete[] w;
      w = nullptr;
    }
  if (cw != nullptr)
    {
      delete[] cw;
      cw = nullptr;
    }
}


// normalize the weights and prepare cumulative weights
void intra_step_distribution::norm_w_prep_cw()
{
  double z = 0.0;

  // determine the normalization constant
  for (int i=0; i<N; i++)
    {
      z += w[i];
    }

  // normalize the weights
  for (int i=0; i<N; i++)
    {
      w[i] = w[i]/z;
    }

  // calculate the cumulative weights
  for (int i=0; i<N; i++)
    {
      if (i == 0)
	{
	  cw[i] = w[i];
	}
      else
	{
	  cw[i] = cw[i-1] + w[i];
	}
    }
}


// print the distribution
void intra_step_distribution::print_distribution()
{
  for (int i=0; i<N; i++)
    {
      std::cout << "i=" << i << "," << w[i] << "," << cw[i] << std::endl;
    }
}


// sample from the distribution
int intra_step_distribution::get_k(double r, int k_max)
{
  int k = 0;
  double rc = r*cw[k_max-1];

  // std::cout << "sampling k, rc=" << rc << std::endl;
  // std::cout << r << "," << k_max << "," << cw[k_max-1];

  for (int i=0; i<k_max; i++)
    {
      if (rc < cw[i])
	{
	  k = i;
	  break;
	}
    }
  return k;
}


// prepare a uniform distribution for future step sampling
void intra_step_distribution::uniform_dist(int k_max)
{
  set_N(k_max);
  uniform_w();
  norm_w_prep_cw();
  std::cout << "uniform step distribution" << std::endl;
  print_distribution();
}


// set the weights according to a uniform distribution
void intra_step_distribution::uniform_w()
{
  for (int i=0; i<N; i++)
    {
      w[i] = 1.0;
    }
}


// prepare a truncated Poisson distribution for future step sampling
void intra_step_distribution::poisson_dist(double l, int k_max)
{
  set_N(k_max);
  poisson_w(l);
  norm_w_prep_cw();
  std::cout << "Poisson step distribution" << std::endl;
  print_distribution();
}


// set the weights according to a Poisson distribution
void intra_step_distribution::poisson_w(double l)
{
  for (int i=0; i<N; i++)
    {
      w[i] = exp(i*log(l)-l-lgamma(i+1));
    }
}
