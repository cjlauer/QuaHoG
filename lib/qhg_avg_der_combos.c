#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include <complex.h>
#include <qhg_idx.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_der_correlator.h>
#include <math.h>
#define VG(v, g) (v*NS*NS+g)
#define DD(mu,nu) (ND*mu+nu)
#define DDD(mu,nu,rho) (ND*ND*mu+ND*nu+rho)

/* 
   Calculates average of 2nd derivative correlators of the form
   O_ijk = g_i D_j D_k over permutations of two indices and 0 
*/

void qhg_avg_der_combos_O_0_i_j(_Complex double *out, _Complex double **in, 
				int i, int j, size_t vol_size) {

  for(size_t v=0; v<vol_size; v++){

    out[v] = ( in[DD(i,j)][VG(v,0+2)] 
	       + in[DD(j,i)][VG(v,0+2)] 
	       + in[DD(0,j)][VG(v,i+2)] 
	       + in[DD(j,0)][VG(v,i+2)] 
	       + in[DD(0,i)][VG(v,j+2)] 
	       + in[DD(i,0)][VG(v,j+2)] ) / 6;
    
  }

  return;
}

/*
  Calculate average of 3rd derivative correlators of the form 
  O_{mu mu nu nu} = g_mu D_mu D_nu D_nu over permutations of two indices
*/
void qhg_avg_der_combos_O_mu_mu_nu_nu(_Complex double *out, _Complex double **in, 
				      int mu, int nu, size_t vol_size) {
  
  for(size_t v=0; v<vol_size; v++){

    out[v] = ( in[DDD(mu,nu,nu)][VG(v,mu+2)] 
	       + in[DDD(nu,mu,nu)][VG(v,mu+2)]
	       + in[DDD(nu,nu,mu)][VG(v,mu+2)]
	       + in[DDD(mu,mu,nu)][VG(v,nu+2)]
	       + in[DDD(mu,nu,mu)][VG(v,nu+2)]
	       + in[DDD(nu,mu,mu)][VG(v,nu+2)] ) / 6;

  }

  return;
}


/*
  Calculate average of 3rd derivative correlators of the form 
  O_{mu nu rho sig} g_mu D_nu D_rho D_sig over permutations 
  of four indices without any pairs of equal indices
*/
void qhg_avg_der_combos_O_mu_nu_rho_sig(_Complex double *out, _Complex double **in, 
					int mu, int nu, int rho, int sig, 
					size_t vol_size) {

  for(size_t v=0; v<vol_size; v++){

    out[v] = ( in[DDD(nu,rho,sig)][VG(v,mu+2)]
	       + in[DDD(nu,sig,rho)][VG(v,mu+2)]
	       + in[DDD(rho,nu,sig)][VG(v,mu+2)]
	       + in[DDD(rho,sig,nu)][VG(v,mu+2)]
	       + in[DDD(sig,nu,rho)][VG(v,mu+2)]
	       + in[DDD(sig,rho,nu)][VG(v,mu+2)]
	       + in[DDD(mu,rho,sig)][VG(v,nu+2)]
	       + in[DDD(mu,sig,rho)][VG(v,nu+2)]
	       + in[DDD(rho,mu,sig)][VG(v,nu+2)]
	       + in[DDD(rho,sig,mu)][VG(v,nu+2)]
	       + in[DDD(sig,mu,rho)][VG(v,nu+2)]
	       + in[DDD(sig,rho,mu)][VG(v,nu+2)]
	       + in[DDD(mu,nu,sig)][VG(v,rho+2)]
	       + in[DDD(mu,sig,nu)][VG(v,rho+2)]
	       + in[DDD(nu,mu,sig)][VG(v,rho+2)]
	       + in[DDD(nu,sig,mu)][VG(v,rho+2)]
	       + in[DDD(sig,mu,nu)][VG(v,rho+2)]
	       + in[DDD(sig,nu,mu)][VG(v,rho+2)]
	       + in[DDD(mu,nu,rho)][VG(v,sig+2)]
	       + in[DDD(mu,rho,nu)][VG(v,sig+2)]
	       + in[DDD(nu,mu,rho)][VG(v,sig+2)]
	       + in[DDD(nu,rho,mu)][VG(v,sig+2)]
	       + in[DDD(rho,mu,nu)][VG(v,sig+2)]
	       + in[DDD(rho,nu,mu)][VG(v,sig+2)] ) / 24;

  }

  return;
}

qhg_der_correlator qhg_avg_der_combos_der2(qhg_der_correlator in) {

  size_t site_size = 1;
  int ncorr = 3;

  qhg_der_correlator out = qhg_averaged_der_correlator_init(site_size, in.lat, in.mom_list, ncorr, in.der_order);

  int c = 0;
  for(int i=1; i<ND; i++)
    for(int j=i+1; j<ND; j++){
      qhg_avg_der_combos_O_0_i_j(out.C[c], in.C, i, j, in.vol_size);
      c++;
    }

  if(in.mom_list == NULL) {
    for(int i=0; i<ND; i++)
      out.origin[i] = in.origin[i];
    out.mom_list = NULL;
  } else {
    out.origin[0] = in.origin[0];
    out.mom_list = in.mom_list;
  }

  out.dt = in.dt;
  out.flav = in.flav;
  out.proj = in.proj;

  return out;
}

qhg_der_correlator qhg_avg_der_combos_der3(qhg_der_correlator in, int *mom_vec) {

  size_t site_size = 1;
  int ncorr = ( mom_vec[0] == 0 || mom_vec[1] == 0 || mom_vec[2] == 0 
		? 6 : 7);

  qhg_der_correlator out = qhg_averaged_der_correlator_init(site_size, 
							    in.lat, 
							    in.mom_list, 
							    ncorr, 
							    in.der_order);

  int c = 0;

  for(int mu=0; mu<ND; mu++)
    for(int nu=mu+1; nu<ND; nu++){
      qhg_avg_der_combos_O_mu_mu_nu_nu(out.C[c], in.C, mu, nu, in.vol_size);
      c++;
    }

  if( mom_vec[0] != 0 && mom_vec[1] != 0 && mom_vec[2] != 0 )
    qhg_avg_der_combos_O_mu_nu_rho_sig(out.C[c], in.C, 
				       0, 1, 2, 3, in.vol_size);

  if(in.mom_list == NULL) {
    for(int i=0; i<ND; i++)
      out.origin[i] = in.origin[i];
    out.mom_list = NULL;
  } else {
    out.origin[0] = in.origin[0];
    out.mom_list = in.mom_list;
  }

  out.dt = in.dt;
  out.flav = in.flav;
  out.proj = in.proj;

  return out;
}

qhg_der_correlator qhg_avg_der_combos(qhg_der_correlator in, int *mom_vec) {

  qhg_der_correlator out;

  if(in.der_order == 2)
    out = qhg_avg_der_combos_der2(in);
  else if(in.der_order == 3)
    out = qhg_avg_der_combos_der3(in, mom_vec);
  else {
    fprintf(stderr, "ERROR: derivative order not supported for averaging.\n");
  }

  return out;
}
