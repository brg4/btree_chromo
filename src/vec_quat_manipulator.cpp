#include <vec_quat_manipulator.hpp>

// constructor
vec_quat_manipulator::vec_quat_manipulator()
{
}

// destructor
vec_quat_manipulator::~vec_quat_manipulator()
{
}

// vector functions

vec vec_quat_manipulator::v_null()
{
  vec v;
  v.x = 0.0;
  v.y = 0.0;
  v.z = 0.0;
  return v;
}

vec vec_quat_manipulator::v_new(double x, double y, double z)
{
  vec v;
  v.x = x;
  v.y = y;
  v.z = z;
  return v;
}

double vec_quat_manipulator::v_dot(vec v, vec w)
{
  return v.x*w.x + v.y*w.y + v.z*w.z;
}

double vec_quat_manipulator::v_L2(vec v)
{
  return v_dot(v,v);
}

vec vec_quat_manipulator::v_norm(vec v)
{
  return v_ax(1.0/sqrt(v_L2(v)),v);
}

vec vec_quat_manipulator::v_inv(vec v)
{
  return v_ax(-1.0,v);
}

vec vec_quat_manipulator::v_ax(double a, vec v)
{
  vec ax;
  ax.x = a*v.x;
  ax.y = a*v.y;
  ax.z = a*v.z;
  return ax;
}

vec vec_quat_manipulator::v_xpy(vec v, vec w)
{
  vec xpy;
  xpy.x = v.x + w.x;
  xpy.y = v.y + w.y;
  xpy.z = v.z + w.z;
  return xpy;
}

vec vec_quat_manipulator::v_axpy(double a, vec v, vec w)
{ 
  return v_xpy(v_ax(a,v),w);
}

vec vec_quat_manipulator::v_cross(vec v, vec w)
{
  vec c;
  c.x = v.y*w.z - v.z*w.y;
  c.y = v.z*w.x - v.x*w.z;
  c.z = v.x*w.y - v.y*w.x;
  return c;
}

vec vec_quat_manipulator::v_linterp(double s, vec v, vec w)
{
  return v_axpy(s,v_xpy(v,v_inv(w)),w);
}

void vec_quat_manipulator::v_array_from_vector(double **&x, int &N, std::vector<vec> &vs)
{

  N = static_cast<int>(vs.size());
  
  x = new double*[N];
  for (int i=0; i<N; i++)
    {
      x[i] = new double[3];
      x[i][0] = vs[i].x;
      x[i][1] = vs[i].y;
      x[i][2] = vs[i].z;
    }
  
}

void vec_quat_manipulator::v_vector_from_array(std::vector<vec> &vs, double **&x, int &N)
{

  vs.clear();
  
  for (int i=0; i<N; i++)
    {
      vs.push_back(v_new(x[i][0],x[i][1],x[i][2]));
    }
  
}

// quaternion functions

quat vec_quat_manipulator::q_null()
{
  quat q;
  q.w = 0.0;
  q.v = v_null();
  return q;
}

double vec_quat_manipulator::q_dot(quat q, quat p)
{
  return q.w*p.w + v_dot(q.v,p.v);
}

double vec_quat_manipulator::q_L2(quat q)
{
  return q_dot(q,q);
}

quat vec_quat_manipulator::q_norm(quat q)
{
  return q_ax(1.0/sqrt(q_L2(q)),q);
}

quat vec_quat_manipulator::q_ax(double a, quat q)
{
  quat ax;
  ax.w = a*q.w;
  ax.v = v_ax(a,q.v);
  return ax;
}

quat vec_quat_manipulator::q_xpy(quat q, quat p)
{
  quat xpy;
  xpy.w = q.w + p.w;
  xpy.v = v_xpy(q.v,p.v);
  return xpy;
}

quat vec_quat_manipulator::q_axpy(double a, quat q, quat p)
{ 
  return q_ax(a,q_xpy(q,p));
}

quat vec_quat_manipulator::q_mult(quat q, quat p)
{
  quat qmp;
  qmp.w = q.w*p.w - v_dot(q.v,p.v);
  qmp.v = v_cross(q.v,p.v);
  qmp.v = v_axpy(p.w,q.v,qmp.v);
  qmp.v = v_axpy(q.w,p.v,qmp.v);
  return qmp;
}

quat vec_quat_manipulator::q_conj(quat q)
{
  quat qc = q;
  qc.v = v_ax(-1.0,qc.v);
  return qc;
}

quat vec_quat_manipulator::q_inv(quat q)
{
  return q_ax(1.0/q_L2(q),q_conj(q));
}


// conversions

// convert a vector to a quaternion
quat vec_quat_manipulator::v_to_q(vec v)
{
  quat q = q_null();
  q.v = v;
  return q;
}

// convert a quaternion to a vector
vec vec_quat_manipulator::q_to_v(quat q)
{
  return q.v;
}
