#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <complex.h>
#include <qhg_defs.h>
#include <qhg_idx.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_xchange_utils.h>

/*
  This will swap every site with it's neighbor in +dir 

  1. If sign == +1, this will begin from the first site in direction
  "dir", and move forward in direction "dir". In the end the first
  slice will be the last
  
  2. If sign == -1, this will begin from the last slice in direction "dir"
  and move backwards in dir. In the end the last slice will be the first.  
 */
static void
shift_locally(_Complex double *c, unsigned long int **nn, size_t site_size, int dims[ND], int dir, int sign)
{
  unsigned long int lvol = gf.lat->lvol;
  unsigned long int **nn = gf.lat->nn;
  _Complex double swap[site_size];  

  unsigned long int v0, v1;
  for(int t=0; t<dims[0]; t++)
    for(int x=0; x<dims[1]; x++)
      for(int y=0; y<dims[2]; y++)
	for(int z=0; z<dims[3]; z++) {
	  unsigned long int xv[] = {t,x,y,z};
	  if(sign == -1) {
            xv[dir] = dims[dir] - xv[dir] - 1;
            v0 = IDX(xv, dims);
            v1 = nn[dir+ND][v0];
	  } else {
            v0 = IDX(xv, dims);
            v1 = nn[dir][v0];
          }
	  memcpy(swap, &c[v0*site_size], sizeof(_Complex double)*site_size);
	  memcpy(&c[v0*site_size], &c[v1*site_size], sizeof(_Complex double)*site_size);
	  memcpy(&c[v1*site_size], swap, sizeof(_Complex double)*site_size);	    
	}
  return;
}

void
qhg_gauge_field_shift(qhg_gauge_field gf, int shifts_in[ND])
{
  qhg_lattice *lat = gf.lat;
  unsigned long int lvol = gf.lat->lvol;
  unsigned long int **nn = gf.lat->nn;
  int *dims = lat->dims;
  int shifts[ND];
  /* First ensure shift is within bounds */
  for(int dir=0; dir<ND; dir++)
    shifts[dir] = (shifts_in[dir] + dims[dir]) % dims[dir];
  
  /* Get the shortest absolute distance */
  for(int dir=0; dir<ND; dir++)
    shifts[dir] = shifts[dir] <= dims[dir]/2 ? shifts[dir] : (shifts[dir] - dims[dir]);
  
  int site_size = NC*NC*ND;
  _Complex double *U = gf.field;
  for(int dir=0; dir<ND; dir++) {
    int nn_dir = shifts[dir] < 0 ? dir : (dir+ND);
    for(int i=0; i<abs(shifts[dir]); i++) {
      shift_locally(U, nn, site_size, dims, dir, sign);
      qhg_xchange_gauge(gf);
    }
  }

  return;
}
