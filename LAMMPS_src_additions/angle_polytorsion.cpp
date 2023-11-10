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
#include "angle_polytorsion.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

#include "atom_vec_ellipsoid.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AnglePolytorsion::AnglePolytorsion(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AnglePolytorsion::~AnglePolytorsion()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(kalign);
    memory->destroy(ktwist);
  }
}

/* ---------------------------------------------------------------------- */

void AnglePolytorsion::compute(int eflag, int vflag)
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

  eangle = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x; // position vector
  double **f = atom->f; // force vector
  double **torque = atom->torque;
  int *ellipsoid = atom->ellipsoid;
  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;


  if (!newton_bond)
    error->all(FLERR,"'newton' flag for bonded interactions must be 'on'");

  for (n = 0; n < nanglelist; n++) {
    i = anglelist[n][0]; 
    iP1 = anglelist[n][1]; 
    iDummy = anglelist[n][2]; 
    type = anglelist[n][3];

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

    tI[0] = x[iP1][0] - x[i][0]; // get vector between i and iP1
    tI[1] = x[iP1][1] - x[i][1];
    tI[2] = x[iP1][2] - x[i][2];
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
    
    torque[i][0] += kalign[type]*uIcrosstI[0] + kHI[0];  
    torque[i][1] += kalign[type]*uIcrosstI[1] + kHI[1]; 
    torque[i][2] += kalign[type]*uIcrosstI[2] + kHI[2];
 
    torque[iP1][0] -= kHI[0];  
    torque[iP1][1] -= kHI[1]; 
    torque[iP1][2] -= kHI[2]; 
    

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
      
    f[i][0] += f1[0];
    f[i][1] += f1[1];
    f[i][2] += f1[2];
    
    f[iP1][0] -= f1[0];
    f[iP1][1] -= f1[1];
    f[iP1][2] -= f1[2];
    
    f3[0] = f3[1] = f3[2]  = 0.0;  // for virial calculation

    if (eflag) eangle = kalign[type] * (1 - uIdottI) + ktwist[type] * (1 - cosalphaiplusgammai);

    if (evflag) // tally energy (virial=0 because force=0)
      ev_tally(i,iP1,iDummy,nlocal,newton_bond,eangle,f1,f3,
               dx,dy,dz,0.0,0.0,0.0);

  }
}

/* ---------------------------------------------------------------------- */

void AnglePolytorsion::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(kalign,n+1,"angle:kalign");
  memory->create(ktwist,n+1,"angle:ktwist");

  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AnglePolytorsion::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nangletypes,ilo,ihi,error);

  double kalign_one = utils::numeric(FLERR,arg[1],false,lmp);
  double ktwist_one = utils::numeric(FLERR,arg[2],false,lmp);

  // convert gamma0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    kalign[i] = kalign_one;
    ktwist[i] = ktwist_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ----------------------------------------------------------------------
   used by SHAKE
------------------------------------------------------------------------- */

double AnglePolytorsion::equilibrium_angle(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AnglePolytorsion::write_restart(FILE *fp)
{
  fwrite(&kalign[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ktwist[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AnglePolytorsion::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&kalign[1],sizeof(double),atom->nangletypes,fp);
    fread(&ktwist[1],sizeof(double),atom->nangletypes,fp);
  }
  MPI_Bcast(&kalign[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ktwist[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   used by ComputeAngleLocal
------------------------------------------------------------------------- */

double AnglePolytorsion::single(int type, int i, int iP1, int iDummy)
{
  double **x = atom->x; // position vector
  int *ellipsoid = atom->ellipsoid;
  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;

  double delx = x[iP1][0] - x[i][0];
  double dely = x[iP1][1] - x[i][1];
  double delz = x[iP1][2] - x[i][2];

  domain->minimum_image(delx,dely,delz);

  double  *quatP1 = bonus[ellipsoid[iP1]].quat;
  double fP1[3],vP1[3],uP1[3];
  fP1[0] = quatP1[0]*quatP1[0] + quatP1[1]*quatP1[1] - quatP1[2]*quatP1[2] - quatP1[3]*quatP1[3];
  fP1[1] = 2.0*( quatP1[1]*quatP1[2] + quatP1[0]*quatP1[3] );
  fP1[2] = 2.0*( quatP1[1]*quatP1[3] - quatP1[0]*quatP1[2] );
  vP1[0] = 2.0*( quatP1[1]*quatP1[2] - quatP1[0]*quatP1[3] );
  vP1[1] = quatP1[0]*quatP1[0] - quatP1[1]*quatP1[1] + quatP1[2]*quatP1[2] - quatP1[3]*quatP1[3];
  vP1[2] = 2.0*( quatP1[2]*quatP1[3] + quatP1[0]*quatP1[1] );
  uP1[0] = 2.0*( quatP1[1]*quatP1[3] + quatP1[0]*quatP1[2] );
  uP1[1] = 2.0*( quatP1[2]*quatP1[3] - quatP1[0]*quatP1[1] );
  uP1[2] = quatP1[0]*quatP1[0] - quatP1[1]*quatP1[1] - quatP1[2]*quatP1[2] + quatP1[3]*quatP1[3];

  double  *quatI = bonus[ellipsoid[i]].quat;
  double fI[3],vI[3],uI[3];
  fI[0] = quatI[0]*quatI[0] + quatI[1]*quatI[1] - quatI[2]*quatI[2] - quatI[3]*quatI[3];
  fI[1] = 2.0*( quatI[1]*quatI[2] + quatI[0]*quatI[3] );
  fI[2] = 2.0*( quatI[1]*quatI[3] - quatI[0]*quatI[2] );
  vI[0] = 2.0*( quatI[1]*quatI[2] - quatI[0]*quatI[3] );
  vI[1] = quatI[0]*quatI[0] - quatI[1]*quatI[1] + quatI[2]*quatI[2] - quatI[3]*quatI[3];
  vI[2] = 2.0*( quatI[2]*quatI[3] + quatI[0]*quatI[1] );
  uI[0] = 2.0*( quatI[1]*quatI[3] + quatI[2]*quatI[0] );
  uI[1] = 2.0*( quatI[2]*quatI[3] - quatI[1]*quatI[0] );
  uI[2] = quatI[0]*quatI[0] - quatI[1]*quatI[1] - quatI[2]*quatI[2] + quatI[3]*quatI[3];


  double inv_uIdotuiP1plus1 = 1.0/(uI[0]*uP1[0] + uI[1]*uP1[1] + uI[2]*uP1[2] + 1);
  double cosalphaiplusgammai = (fP1[0]*fI[0] + fP1[1]*fI[1] + fP1[2]*fI[2] + 
				vP1[0]*vI[0] + vP1[1]*vI[1] + vP1[2]*vI[2]) * inv_uIdotuiP1plus1;

  double bI2 = delx*delx + dely*dely + delz*delz;
  double cosphii = uI[0]*delx + uI[1]*dely + uI[2]*delz;
  cosphii/=sqrt(bI2);

  return kalign[type]*(1-cosphii) +  ktwist[type]*(1-cosalphaiplusgammai); // energy
}
