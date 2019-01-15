#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>
#include <complex.h>
#include <qhg_defs.h>
#include <qhg_su3_ops.h>
#include <qhg_xchange_gauge.h>
#include <qhg_alloc.h>

#define G(v,mu) (((v)*ND + mu)*NC*NC)

qhg_gluon_loop
qhg_gluon_loop_init(qhg_lattice *lat)
{
  unsigned long int lvol = lat->lvol;

  qhg_gluon_loop gl;
  gl.loop = qhg_alloc(lvol*sizeof(_Complex double));
    
  gl.lat = lat;
  return gl;
}

void
qhg_gluon_loop_finalize(qhg_gluon_loop gl)
{
  free(gl.loop);
  gl.lat = NULL;
  return;
}

void
qhg_gluon_loop_copy(qhg_gluon_loop y, qhg_gluon_loop x)
{
  y.lat = x.lat;
  memcpy(y.loop, x.loop, x.lat->lvol*sizeof(_Complex double));
  return;
}

static _Complex double
gluon_plaq(_Complex double *U, unsigned long int **nn, unsigned long int v00, int mu, int nu)
{
  double plaq = 0.0;
  _Complex double x[NC*NC];
  _Complex double y[NC*NC];	
  _Complex double *u0, *u1, *u2, *u3;
  unsigned long int vmu = nn[mu][v00];
  unsigned long int vnu = nn[nu][v00];	
  u0 = &U[G(v00, mu)];
  u1 = &U[G(vmu, nu)];
  u2 = &U[G(vnu, mu)];
  u3 = &U[G(v00, nu)];
  su3_mul_uu(x, u0, u1);
  su3_mul_ud(y, x, u2);
  su3_mul_ud(x, y, u3);	
      
  return su3_linalg_trace_u(x);
}

void
qhg_calculate_gluon_loop(qhg_gluon_loop out, qhg_gauge_field gf)
{
  qhg_lattice *lat = gf.lat;
  qhg_comms *comms = lat->comms;

  _Complex double *U = gf.field;

  qhg_xchange_gauge(gf);

  unsigned long int lvol = lat->lvol;
  unsigned long int **nn = lat->nn;

  qhg_gluon_loop gl = qhg_gluon_loop_init(lat);

  for(unsigned long int v00=0; v00<lvol; v00++)
    for(unsigned long int i=1; i<ND; i++){

      gl.loop[v00] += gluon_plaq(U, nn, v00, i, 0);

      for(unsigned long int j=1; j<i; j++)
	gl.loop[v00] -= gluon_plaq(U, nn, v00, j, i);
    
    }

  qhg_gluon_loop_copy(out, gl);

  qhg_gluon_loop_finalize(gl);

  return;
}
  
