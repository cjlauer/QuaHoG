#include <stdlib.h>
#include <string.h>

#include <string.h>
#include <complex.h>
#include <qhg_defs.h>
#include <qhg_idx.h>
#include <qhg_types.h>

void
qhg_point_spinor_field(qhg_spinor_field sp, int coords[ND], int spin, int color)
{
  int *ldims = sp.lat->ldims;
  int *procs = sp.lat->comms->proc_dims;
  unsigned long int lvol = sp.lat->lvol;
  int proc_id = sp.lat->comms->proc_id;
  
  /* Coordinates of proc which holds coords[], in t,x,y,z */
  int c_proc[ND];
  for(int i=0; i<ND; i++)
    c_proc[i] = coords[i]/ldims[i];

  /* Which proc has c_proc[] coordinates in the process grid */
  int s_proc = IDX(c_proc, procs);

  /* The local coords of the point source */
  int lc[ND];
  for(int i=0; i<ND; i++)
    lc[i] = coords[i] % ldims[i];

  unsigned long int lidx = IDX(lc, ldims);

  /* Set source to zero on all procs */
  memset(sp.field, 0, lvol*NC*NS*sizeof(_Complex double));

  /* Set spinor to 1 at [coords,spin,color] */
  if(proc_id == s_proc)
    sp.field[lidx*NC*NS + CS(spin,color)] = 1.0;

  return;
}
