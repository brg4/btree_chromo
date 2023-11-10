/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Chris A Brackley (U Edinburgh) 
   Based on the simulation scheme described in 
       G Chirico and J Langowski, Biopolymers 34 p415-433 (1994)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "omp_compat.h"
#include "angle_polytorsion_omp.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

#include "atom_vec_ellipsoid.h"

#include "suffix.h"
using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AnglePolytorsionOMP::AnglePolytorsionOMP(class LAMMPS *lmp)
  : AnglePolytorsion(lmp), ThrOMP(lmp,THR_ANGLE)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void AnglePolytorsionOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = neighbor->nanglelist;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, cvatom, thr);

    if (inum > 0) {
      if (evflag) {
        if (eflag) {
          if (force->newton_bond) eval<1,1,1>(ifrom, ito, thr);
          else eval<1,1,0>(ifrom, ito, thr);
        } else {
          if (force->newton_bond) eval<1,0,1>(ifrom, ito, thr);
          else eval<1,0,0>(ifrom, ito, thr);
        }
      } else {
        if (force->newton_bond) eval<0,0,1>(ifrom, ito, thr);
        else eval<0,0,0>(ifrom, ito, thr);
      }
    }
    thr->timer(Timer::BOND);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int NEWTON_BOND>
void AnglePolytorsionOMP::eval(int nfrom, int nto, ThrData * const thr)
{
  int i,iP1,iDummy,n,type;
  double eangle,f1[3],f3[3],dx,dy,dz;
  double uI[3],uP1[3],  // axis of the atoms
    fI[3],fP1[3],
    vI[3],vP1[3];

  double *quatI, *quatP1;

  double uIdotuiP1plus1, inv_uIdotuiP1plus1, cosalphaiplusgammai;
  double fIcrossfP1[3], vIcrossvP1[3], uIcrossuP1[3], HI[3], kHI[3];
  double tI[3], inv_bI, bI2;
  double uIcrosstI[3];

  double GI[3],uIdottI;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  dbl3_t * _noalias const torque = (dbl3_t *) thr->get_torque()[0];
  // int *ellipsoid = atom->ellipsoid;
  // AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  const int * _noalias const ellipsoid = atom->ellipsoid;
  const int4_t * _noalias const anglelist = (int4_t *) neighbor->anglelist[0];
  const int nlocal = atom->nlocal;

  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  eangle = 0.0;

  if (!NEWTON_BOND)
    error->all(FLERR,"'newton' flag for bonded interactions must be 'on'");

  for (n = nfrom; n < nto; n++) {
    i = anglelist[n].a; 
    iP1 = anglelist[n].b; 
    iDummy = anglelist[n].c; 
    type = anglelist[n].t;

    /* ----------------------------------------------
       Get u,f,v axes of beads i and iP1
       ---------------------------------------------- */
    quatP1 = bonus[ellipsoid[iP1]].quat;
    fP1[0] = quatP1[0]*quatP1[0] + quatP1[1]*quatP1[1] - quatP1[2]*quatP1[2] - quatP1[3]*quatP1[3];
    fP1[1] = 2.0*( quatP1[1]*quatP1[2] + quatP1[0]*quatP1[3] );
    fP1[2] = 2.0*( quatP1[1]*quatP1[3] - quatP1[0]*quatP1[2] );
    vP1[0] = 2.0*( quatP1[1]*quatP1[2] - quatP1[0]*quatP1[3] );
    vP1[1] = quatP1[0]*quatP1[0] - quatP1[1]*quatP1[1] + quatP1[2]*quatP1[2] - quatP1[3]*quatP1[3];
    vP1[2] = 2.0*( quatP1[2]*quatP1[3] + quatP1[0]*quatP1[1] );
    uP1[0] = 2.0*( quatP1[1]*quatP1[3] + quatP1[0]*quatP1[2] );
    uP1[1] = 2.0*( quatP1[2]*quatP1[3] - quatP1[0]*quatP1[1] );
    uP1[2] = quatP1[0]*quatP1[0] - quatP1[1]*quatP1[1] - quatP1[2]*quatP1[2] + quatP1[3]*quatP1[3];

    quatI = bonus[ellipsoid[i]].quat;
    fI[0] = quatI[0]*quatI[0] + quatI[1]*quatI[1] - quatI[2]*quatI[2] - quatI[3]*quatI[3];
    fI[1] = 2.0*( quatI[1]*quatI[2] + quatI[0]*quatI[3] );
    fI[2] = 2.0*( quatI[1]*quatI[3] - quatI[0]*quatI[2] );
    vI[0] = 2.0*( quatI[1]*quatI[2] - quatI[0]*quatI[3] );
    vI[1] = quatI[0]*quatI[0] - quatI[1]*quatI[1] + quatI[2]*quatI[2] - quatI[3]*quatI[3];
    vI[2] = 2.0*( quatI[2]*quatI[3] + quatI[0]*quatI[1] );
    uI[0] = 2.0*( quatI[1]*quatI[3] + quatI[0]*quatI[2] );
    uI[1] = 2.0*( quatI[2]*quatI[3] - quatI[0]*quatI[1] );
    uI[2] = quatI[0]*quatI[0] - quatI[1]*quatI[1] - quatI[2]*quatI[2] + quatI[3]*quatI[3];

    /* ----------------------------------------------
       Twist torque
       ---------------------------------------------- */

    uIdotuiP1plus1 = uI[0]*uP1[0] + uI[1]*uP1[1] + uI[2]*uP1[2] + 1;
    inv_uIdotuiP1plus1 = 1.0/uIdotuiP1plus1;
    cosalphaiplusgammai = (fP1[0]*fI[0] + fP1[1]*fI[1] + fP1[2]*fI[2] + 
			   vP1[0]*vI[0] + vP1[1]*vI[1] + vP1[2]*vI[2]) * inv_uIdotuiP1plus1;

    fIcrossfP1[0] = fI[1]*fP1[2] - fI[2]*fP1[1];
    fIcrossfP1[1] = fI[2]*fP1[0] - fI[0]*fP1[2];
    fIcrossfP1[2] = fI[0]*fP1[1] - fI[1]*fP1[0];
    vIcrossvP1[0] = vI[1]*vP1[2] - vI[2]*vP1[1];
    vIcrossvP1[1] = vI[2]*vP1[0] - vI[0]*vP1[2];
    vIcrossvP1[2] = vI[0]*vP1[1] - vI[1]*vP1[0];
    uIcrossuP1[0] = uI[1]*uP1[2] - uI[2]*uP1[1];
    uIcrossuP1[1] = uI[2]*uP1[0] - uI[0]*uP1[2];
    uIcrossuP1[2] = uI[0]*uP1[1] - uI[1]*uP1[0];

    HI[0] = inv_uIdotuiP1plus1 * 
      (  fIcrossfP1[0] +  vIcrossvP1[0] - cosalphaiplusgammai*uIcrossuP1[0]);
    HI[1] = inv_uIdotuiP1plus1 * 
      (  fIcrossfP1[1] +  vIcrossvP1[1] - cosalphaiplusgammai*uIcrossuP1[1]);
    HI[2] = inv_uIdotuiP1plus1 * 
      (  fIcrossfP1[2] +  vIcrossvP1[2] - cosalphaiplusgammai*uIcrossuP1[2]);

    kHI[0] = ktwist[type] * HI[0];
    kHI[1] = ktwist[type] * HI[1];
    kHI[2] = ktwist[type] * HI[2];

    /* ----------------------------------------------
       Allignment torque
       ---------------------------------------------- */

    tI[0] = x[iP1].x - x[i].x; // get vector between i and iP1
    tI[1] = x[iP1].y - x[i].y;
    tI[2] = x[iP1].z - x[i].z;
    dx=tI[0]; // store bead separations for virial calculation later
    dy=tI[1];
    dz=tI[2];
    bI2 = tI[0]*tI[0] + tI[1]*tI[1] + tI[2]*tI[2];
    inv_bI = 1.0/sqrt(bI2);
    tI[0]*=inv_bI;  // make tI a unit vector
    tI[1]*=inv_bI;
    tI[2]*=inv_bI;

    uIcrosstI[0] = uI[1]*tI[2] - uI[2]*tI[1];
    uIcrosstI[1] = uI[2]*tI[0] - uI[0]*tI[2];
    uIcrosstI[2] = uI[0]*tI[1] - uI[1]*tI[0];

    /* ----------------------------------------------
       Total torque
       ---------------------------------------------- */
    
    torque[i].x += kalign[type]*uIcrosstI[0] + kHI[0];  
    torque[i].y += kalign[type]*uIcrosstI[1] + kHI[1]; 
    torque[i].z += kalign[type]*uIcrosstI[2] + kHI[2];
 
    torque[iP1].x -= kHI[0];  
    torque[iP1].y -= kHI[1]; 
    torque[iP1].z -= kHI[2]; 
    

    /* ----------------------------------------------
       Force
       ---------------------------------------------- */

    uIdottI = uI[0]*tI[0] + uI[1]*tI[1] + uI[2]*tI[2];

    GI[0] = inv_bI * ( uIdottI*tI[0] - uI[0] );
    GI[1] = inv_bI * ( uIdottI*tI[1] - uI[1] );
    GI[2] = inv_bI * ( uIdottI*tI[2] - uI[2] );
   
    f1[0] = kalign[type] * GI[0] ; 
    f1[1] = kalign[type] * GI[1] ;
    f1[2] = kalign[type] * GI[2] ;

    // apply force to each of 3 atoms

    if (NEWTON_BOND || i < nlocal) {
      f[i].x += f1[0];
      f[i].y += f1[1];
      f[i].z += f1[2];
    }

    if (NEWTON_BOND || iP1 < nlocal) {
      f[iP1].x -= f1[0];
      f[iP1].y -= f1[1];
      f[iP1].z -= f1[2];
    }

    f3[0] = f3[1] = f3[2]  = 0.0;  // for virial calculation

    if (EFLAG) eangle = kalign[type] * (1 - uIdottI) + ktwist[type] * (1 - cosalphaiplusgammai);

    if (EVFLAG) // tally energy (virial=0 because force=0)
      ev_tally_thr(this,i,iP1,iDummy,nlocal,NEWTON_BOND,eangle,f1,f3,
		   dx,dy,dz,0.0,0.0,0.0,thr);
  }
}
