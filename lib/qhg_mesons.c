#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>
#include <qhg_idx.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_correlator.h>
#include <qhg_prop_ops.h>
#include <qhg_prop_gammas.h>
#include <qhg_meson_defs.h>

qhg_correlator
qhg_mesons(qhg_spinor_field sp_u[NS*NC], qhg_spinor_field sp_d[NS*NC], int source_coords[ND])
{
  qhg_lattice *lat = sp_u[0].lat;
  qhg_correlator corr = qhg_correlator_init(NCHANNELS, lat);
  for(int i=0; i<ND; i++)
    corr.origin[i] = source_coords[i];
  int tsrc = corr.origin[0];  
  unsigned long int lvol = lat->lvol;
  unsigned long int lv3 = lat->lv3;
  int lt = lat->ldims[0];
  int t0 = lat->ldims[0]*lat->comms->proc_coords[0];  

#ifdef QHG_OMP
#pragma omp parallel for
#endif
  for(unsigned long int v=0; v<lvol; v++) {
    _Complex double U[NS*NC][NS*NC];
    _Complex double D[NS*NC][NS*NC];    
    _Complex double C[NS*NC][NS*NC];    
    _Complex double W[NS*NC][NS*NC];
    _Complex double V[NS*NC][NS*NC];        

    prop_load(U, sp_u, v);
    prop_load(D, sp_d, v);

    int t = (int)(v/lv3);
    int gt = t + t0;
    if(gt < tsrc) {
      
      switch(sp_u[0].bc) {
      case PERIODIC:
	break;
      case ANTIPERIODIC:
	prop_scale(-1, U);
	break;
      }
    
      switch(sp_d[0].bc) {
      case PERIODIC:
	break;
      case ANTIPERIODIC:
	prop_scale(-1, U);
	break;
      }
      
    }
    
    for(int igamma=0; igamma<NGAMMAS; igamma++) {
      _Complex double (*P[2])[NC*NS] = {U, D};
      for(int iflav=0; iflav<NFLAV; iflav++) {
	/*
	 * Sign flips are for consistency with libqcd
	 */
	switch(igamma) {
	case 0: /* 1 */
	  prop_mul_gd(C, P[iflav], P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = prop_trace(C);
	  break;
	case 1: /* g5 */
	  prop_g5_G(W, P[iflav]);
	  prop_G_g5(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = -prop_trace(C);
	  break;
	case 2: /* gx */
	  prop_gx_G(W, P[iflav]);
	  prop_G_gx(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = -prop_trace(C);
	  break;
	case 3: /* gy */
	  prop_gy_G(W, P[iflav]);
	  prop_G_gy(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = -prop_trace(C);
	  break;
	case 4: /* gz */
	  prop_gz_G(W, P[iflav]);
	  prop_G_gz(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = -prop_trace(C);
	  break;
	case 5: /* gt */
	  prop_gt_G(W, P[iflav]);
	  prop_G_gt(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = prop_trace(C);
	  break;
	case 6: /* g5gx */
	  prop_g5gx_G(W, P[iflav]);
	  prop_G_g5gx(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = -prop_trace(C);
	  break;
	case 7: /* g5gy */
	  prop_g5gy_G(W, P[iflav]);
	  prop_G_g5gy(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = -prop_trace(C);
	  break;
	case 8: /* g5gz */
	  prop_g5gz_G(W, P[iflav]);
	  prop_G_g5gz(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = -prop_trace(C);
	  break;
	case 9: /* g5gt */
	  prop_g5gt_G(W, P[iflav]);
	  prop_G_g5gt(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = prop_trace(C);
	  break;	  
	}
      }
    }
  }
  
  corr.mom_list = NULL;
  return corr;
}
  

qhg_correlator
qhg_mesons_pseudoscalar(qhg_spinor_field sp1[NS*NC], qhg_spinor_field sp2[NS*NC], int source_coords[ND])
{
  qhg_lattice *lat = sp1[0].lat;
  qhg_correlator corr = qhg_correlator_init(1, lat);
  for(int i=0; i<ND; i++)
    corr.origin[i] = source_coords[i];
  int tsrc = corr.origin[0];  
  unsigned long int lvol = lat->lvol;
  unsigned long int lv3 = lat->lv3;
  int lt = lat->ldims[0];
  int t0 = lat->ldims[0]*lat->comms->proc_coords[0];  

#ifdef QHG_OMP
#pragma omp parallel for
#endif
  for(unsigned long int v=0; v<lvol; v++) {
    _Complex double P1[NS*NC][NS*NC];
    _Complex double P2[NS*NC][NS*NC];
    _Complex double C[NS*NC][NS*NC];    

    prop_load(P1, sp1, v);
    prop_load(P2, sp2, v);

    int t = (int)(v/lv3);
    int gt = t + t0;
    if(gt < tsrc) {
      
      switch(sp1[0].bc) {
      case PERIODIC:
	break;
      case ANTIPERIODIC:
	prop_scale(-1, P1);
	break;
      }
    
      switch(sp2[0].bc) {
      case PERIODIC:
	break;
      case ANTIPERIODIC:
	prop_scale(-1, P2);
	break;
      }
      
    }
    
    prop_mul_gd(C, P1, P2);
    corr.C[v] = prop_trace(C);
  }
  
  corr.mom_list = NULL;
  return corr;
}
  
