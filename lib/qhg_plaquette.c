#include <stdlib.h>

#include <mpi.h>
#include <complex.h>
#include <qhg_defs.h>
#include <qhg_su3_ops.h>
#include <qhg_xchange_gauge.h>

#define G(v,mu) (((v)*ND + mu)*NC*NC)

static double
plaquette_00(_Complex double *U, unsigned long int **nn, unsigned long int v00)
{
  double plaq = 0.0;
  _Complex double x[NC*NC];
  _Complex double y[NC*NC];	
  _Complex double *u0, *u1, *u2, *u3;
  for(int mu=0; mu<ND; mu++) 
    for(int nu=mu+1; nu<ND; nu++) {
      unsigned long int vmu = nn[mu][v00];
      unsigned long int vnu = nn[nu][v00];	
      u0 = &U[G(v00, mu)];
      u1 = &U[G(vmu, nu)];
      u2 = &U[G(vnu, mu)];
      u3 = &U[G(v00, nu)];
      su3_mul_uu(x, u0, u1);
      su3_mul_ud(y, x, u2);
      su3_mul_ud(x, y, u3);	
      
      plaq += creal(su3_linalg_trace_u(x));
    }
  return plaq;
}

static double
plaquette_01(_Complex double *U, unsigned long int **nn, unsigned long int v00)
{
  double plaq = 0.0;
  _Complex double x[NC*NC];
  _Complex double y[NC*NC];	
  _Complex double *u0, *u1, *u2, *u3;
  for(int mu=0; mu<ND; mu++) 
    for(int nu=mu+1; nu<ND; nu++) {
      unsigned long int vp0 = nn[mu][v00];
      unsigned long int vpm = nn[mu][nn[nu+ND][v00]];
      unsigned long int v0m = nn[nu+ND][v00];
      u0 = &U[G(v00, mu)];
      u1 = &U[G(vpm, nu)];
      u2 = &U[G(v0m, mu)];
      u3 = &U[G(v0m, nu)];
      su3_mul_ud(x, u0, u1);
      su3_mul_ud(y, x, u2);
      su3_mul_uu(x, y, u3);	
	
      plaq += creal(su3_linalg_trace_u(x));
    }
  return plaq;
}

static double
plaquette_10(_Complex double *U, unsigned long int **nn, unsigned long int v00)
{
  double plaq = 0.0;
  _Complex double x[NC*NC];
  _Complex double y[NC*NC];	
  _Complex double *u0, *u1, *u2, *u3;
  for(int mu=0; mu<ND; mu++) 
    for(int nu=mu+1; nu<ND; nu++) {
      unsigned long int vm0 = nn[mu+ND][v00];
      unsigned long int vmp = nn[mu+ND][nn[nu][v00]];
      unsigned long int v0p = nn[nu][v00];
      u0 = &U[G(vm0, mu)];
      u1 = &U[G(vm0, nu)];
      u2 = &U[G(vmp, mu)];
      u3 = &U[G(v00, nu)];
      su3_mul_du(x, u0, u1);
      su3_mul_uu(y, x, u2);
      su3_mul_ud(x, y, u3);	
	
      plaq += creal(su3_linalg_trace_u(x));
    }
  return plaq;
}

static double
plaquette_11(_Complex double *U, unsigned long int **nn, unsigned long int v00)
{
  double plaq = 0;
  _Complex double x[NC*NC];
  _Complex double y[NC*NC];	
  _Complex double *u0, *u1, *u2, *u3;
  for(int mu=0; mu<ND; mu++) 
    for(int nu=mu+1; nu<ND; nu++) {
      unsigned long int vm0 = nn[mu+ND][v00];
      unsigned long int vmm = nn[mu+ND][nn[nu+ND][v00]];
      unsigned long int v0m = nn[nu+ND][v00];
      u0 = &U[G(vm0, mu)];
      u1 = &U[G(vmm, nu)];
      u2 = &U[G(vmm, mu)];
      u3 = &U[G(v0m, nu)];
      su3_mul_uu(x, u1, u0);
      su3_mul_du(y, u2, x);
      su3_mul_du(x, u3, y);	

      plaq += creal(su3_linalg_trace_u(x));
    }  
  return plaq;
}

double
qhg_plaquette(qhg_gauge_field gf)
{
  unsigned long int lvol = gf.lat->lvol;
  unsigned long int **nn = gf.lat->nn;

  qhg_xchange_gauge(gf);

  _Complex double *U = gf.field;
  
  double plaq = 0.0;
  for(unsigned long int v00=0; v00<lvol; v00++)
    plaq += plaquette_00(U, nn, v00);

  MPI_Comm comm = gf.lat->comms->comm;

  MPI_Allreduce(MPI_IN_PLACE, &plaq, 1, MPI_DOUBLE, MPI_SUM, comm);

  unsigned long int vol = gf.lat->vol;
  return plaq/(ND*(ND-1)/2)/(double)vol/NC;
}
  
