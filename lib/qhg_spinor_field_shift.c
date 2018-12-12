#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <complex.h>
#include <qhg_defs.h>
#include <qhg_idx.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_xchange_utils.h>
#include <qhg_xchange_spinor.h>

/*
  This will swap every site with it's neighbor in +dir 

  1. If sign == +1, this will begin from the first site in direction
  "dir", and move forward in direction "dir". In the end the first
  slice will be the last
  
  2. If sign == -1, this will begin from the last slice in direction "dir"
  and move backwards in dir. In the end the last slice will be the first.  
 */
static void
shift_locally(_Complex double *c, size_t site_size, qhg_lattice *lat, int dir, int sign)
{
  unsigned long int lvol = lat->lvol;
  unsigned long int **nn = lat->nn;

  if(sign>0) {
    for(unsigned long int v = 0; v<lvol; v++) {
      unsigned long int v0 = v;
      unsigned long int v1 = nn[dir][v0];
      memcpy(&c[v0*site_size], &c[v1*site_size], sizeof(_Complex double)*site_size);
    }
  } else {
    for(unsigned long int v = 0; v<lvol; v++) {
      unsigned long int v0 = (lvol-1-v);
      unsigned long int v1 = nn[dir + ND][v0];
      memcpy(&c[v0*site_size], &c[v1*site_size], sizeof(_Complex double)*site_size);
    }
  }
  return;
}

void
qhg_spinor_field_shift(qhg_spinor_field prop, int shifts_in[ND])
{
  qhg_lattice *lat = prop.lat;
  int *dims = lat->dims;
  unsigned long int **nn = prop.lat->nn;
  int shifts[ND];
  /* First ensure shift is within bounds */
  for(int dir=0; dir<ND; dir++)
    shifts[dir] = (shifts_in[dir] + dims[dir]) % dims[dir];
  
  /* Get the shortest absolute distance */
  for(int dir=0; dir<ND; dir++)
    shifts[dir] = shifts[dir] <= dims[dir]/2 ? shifts[dir] : (shifts[dir] - dims[dir]);
  
  qhg_xchange_spinor(prop);
  size_t site_size = NC*NS;
  for(int dir=0; dir<ND; dir++) {
    for(int i=0; i<abs(shifts[dir]); i++) {
      shift_locally(prop.field, site_size, lat, dir, shifts[dir] < 0 ? -1 : 1); 
      qhg_xchange_spinor(prop);
    }
  }

  return;
}
