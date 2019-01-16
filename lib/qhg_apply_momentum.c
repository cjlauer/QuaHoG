#include <mpi.h>
#include <complex.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_idx.h>
#include <qhg_alloc.h>
#include <qhg_spinor_field.h>

void qhg_apply_momentum(qhg_spinor_field out, qhg_spinor_field in, int kx, int ky, int kz, int sign){

  qhg_lattice *lat = in.lat;

  int *ldims = lat->ldims;
  int lt = ldims[0];
  int lx = ldims[1];
  int ly = ldims[2];
  int lz = ldims[3];  

  int Lx = lat->dims[1];
  int Ly = lat->dims[2];
  int Lz = lat->dims[3];  

  int x0 = in.origin[1];
  int y0 = in.origin[2];
  int z0 = in.origin[3];  

  int *pdims = lat->comms->proc_dims;
  unsigned long int lvol = lat->lvol;

  _Complex double *aux = qhg_alloc(lvol*sizeof(_Complex double));
  memset(aux, '\0', lvol*sizeof(_Complex double));
  for(int x=0; x<lx; x++)
    for(int y=0; y<ly; y++)
      for(int z=0; z<lz; z++) {
	int xx = pcoords[1]*lx + x;
	int yy = pcoords[2]*ly + y;
	int zz = pcoords[3]*lz + z;
	double ph =
	  (double)(Lx + xx - x0)*kx/Lx +
	  (double)(Ly + yy - y0)*ky/Ly +
	  (double)(Lz + zz - z0)*kz/Lz;
	ph = ph*2*M_PI;
	_Complex double ex = cos(ph) + sign * _Complex_I * sin(ph);
	for(int t=0; t<lt; t++) {
	  unsigned long int co[] = {t,x,y,z};
	  unsigned long int v = IDX(co, ldims);
	  aux[v] *= ex;
	}
      }
  qhg_spinor_field_copy(out, aux);
  
  qhg_spinor_field_finalize(aux);

  return;

}
