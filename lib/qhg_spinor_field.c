#include <string.h>

#include <complex.h>
#include <qhg_alloc.h>
#include <qhg_defs.h>
#include <qhg_types.h>

#define D(sign, dir) ((sign)*ND + (dir))

qhg_spinor_field
qhg_spinor_field_init(qhg_lattice *lat, enum qhg_fermion_bc_time bc)
{
  unsigned long int lvol = lat->lvol;
  unsigned long int *bvol = lat->bvol;
#ifdef SPINOR_EDGES
  unsigned long int (*evol)[ND] = lat->evol;  
#endif
  unsigned long int ext_vol = lvol;
  for(int i=0; i<ND; i++) {
    ext_vol += 2*bvol[i];
#ifdef SPINOR_EDGES
    for(int j=i+1; j<ND; j++)
      ext_vol += 4*evol[i][j];
#endif    
  }
  
  qhg_spinor_field sp;
  sp.field = qhg_alloc(ext_vol*NS*NC*sizeof(_Complex double));


  unsigned long int v_offset = lvol;
  for(int d0=0; d0<ND; d0++)
    for(int s0=0; s0<2; s0++)    
      if(bvol[d0] != 0) {
	sp.bnd[D(s0, d0)] = &sp.field[v_offset*NC*NS];
	v_offset += bvol[d0];
      } else {
	sp.bnd[D(s0, d0)] = NULL;	
      }

#ifdef SPINOR_EDGES
  for(int d0=0; d0<ND; d0++)
    for(int s0=0; s0<2; s0++)
      for(int d1=d0+1; d1<ND; d1++)
	if(evol[d0][d1] != 0) {
	  for(int s1=0; s1<2; s1++) {
	    sp.edge[D(s0, d0)][D(s1, d1)] = &sp.field[v_offset*NC*NS];
	    v_offset += evol[d0][d1];
	  }
	} else {
	  for(int s1=0; s1<2; s1++)	      
	    sp.edge[D(s0, d0)][D(s1, d1)] = NULL;	
	}
#endif
    
  sp.lat = lat;
  sp.bc = bc;
  return sp;
}

void
qhg_spinor_field_finalize(qhg_spinor_field sp)
{
  free(sp.field);
  for(int dir=0; dir<2*ND; dir++)
    sp.bnd[dir] = NULL;
  sp.lat = NULL;
  return;
}

void
qhg_spinor_field_copy(qhg_spinor_field y, qhg_spinor_field x)
{
  y.lat = x.lat;
  y.bc = x.bc;
  memcpy(y.field, x.field, x.lat->lvol*NC*NS*sizeof(_Complex double));
  return;
}

#define SQUARE(a) (a)*(a)
#define NORM(a) SQUARE(creal(a))+SQUARE(cimag(a))

double
qhg_spinor_field_normsq(qhg_spinor_field y)
{
  unsigned long int lvol = y.lat->lvol;

  double norm = 0;

  for(unsigned long int v=0; v<lvol; v++)
    for(int sc=0; sc<NS*NC; sc++)
      norm += NORM(y.field[v*NS*NC+sc]);
 
  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, y.lat->comms->comm);
  return norm;
}

double
qhg_spinor_field_axpy(double a, qhg_spinor_field x, qhg_spinor_field y)
{
  unsigned long int lvol = y.lat->lvol;

  double norm = 0;

  for(unsigned long int v=0; v<lvol; v++)
    for(int sc=0; sc<NS*NC; sc++)
      x.field[v*NS*NC+sc] = a*x.field[v*NS*NC+sc] + y.field[v*NS*NC+sc];
}
