#ifndef INCLUDE_VEC_QUAT_MANIPULATOR
#define INCLUDE_VEC_QUAT_MANIPULATOR

#include <cmath>
#include <vector>

struct vec
{
  double x, y, z;
};

struct quat
{
  double w;
  vec v;
};

class vec_quat_manipulator
{
public:

  // constructor and destructor
  vec_quat_manipulator();
  ~vec_quat_manipulator();

  vec v_null();
  vec v_new(double x, double y, double z);
  
  double v_dot(vec v, vec w);
  double v_L2(vec v);
  vec v_norm(vec v);
  vec v_inv(vec v);
  vec v_ax(double a, vec v);
  vec v_xpy(vec v, vec w);
  vec v_axpy(double a, vec v, vec w);
  vec v_cross(vec v, vec w);
  vec v_linterp(double s, vec v, vec w);

  void v_array_from_vector(double **&x, int &N, std::vector<vec> &vs);
  void v_vector_from_array(std::vector<vec> &vs, double **&x, int &N);

  quat q_null();

  double q_dot(quat q, quat p);
  double q_L2(quat q);
  quat q_norm(quat q);
  quat q_inv(quat q);
  quat q_ax(double a, quat q);
  quat q_xpy(quat q, quat p);
  quat q_axpy(double a, quat q, quat p);
  quat q_mult(quat q, quat p);
  quat q_conj(quat q);

  quat v_to_q(vec v);
  vec q_to_v(quat q);

private:
  
};

#endif
