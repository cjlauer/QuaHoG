#include <qhg_types.h>
#include <qhg_idx.h>
#include <qhg_prop_ops.h>
#include <qhg_prop_proj.h>
#include <qhg_prop_gammas.h>
#include <math.h>

void
qhg_phase(qhg_spinor_field out[], qhg_spinor_field in[], int mom[3], int source_coords[4], double sign)
{
  qhg_lattice *lat = out[0].lat;
  int *proc_coords = lat->comms->proc_coords;
  int *dims = lat->dims;
  int *ldims = lat->ldims;
  unsigned long int lvol = lat->lvol;

  double shift[3];
  for(int d=0; d<3; d++)
    shift[d] = dims[d+1] + ldims[d+1]*proc_coords[d+1] - source_coords[d+1];

  for(unsigned long int v=0; v<lvol; v++) {
    int co[4] = CO(v,ldims);
    _Complex double coeff = 0;
    double ph = 0;
    for(int d=0; d<3; d++)
      ph += (mom[d]*(shift[d]+co[d+1]))/dims[d+1];

    ph *= sign*2*M_PI;

    coeff = cos(ph) + _Complex_I * sin(ph);

    for(int cs0=0; cs0<NS*NC; cs0++)
      for(int cs1=0; cs1<NS*NC; cs1++) {
        out[cs0].field[VSC(v, cs1)] = in[cs0].field[VSC(v, cs1)]*coeff;
      }
  }
}
