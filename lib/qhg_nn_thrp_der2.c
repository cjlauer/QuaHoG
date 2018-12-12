#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>
#include <qhg_idx.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_correlator.h>
#include <qhg_correlator_shift.h>
#include <qhg_xchange_spinor.h>
#include <qhg_xchange_gauge.h>
#include <qhg_prop_gammas.h>
#include <qhg_prop_ops.h>
#include <qhg_su3_ops.h>
#include <qhg_io_utils.h>
#include <qhg_nn_thrp_der2_defs.h>

static void
prop_bc(enum qhg_fermion_bc_time bc, _Complex double (*p)[NC*NS])
{
  switch(bc) {
  case PERIODIC:
    break;
  case ANTIPERIODIC:
    prop_scale(-1, p);
    break;
  }
  return;
}

qhg_correlator
qhg_nn_thrp_der2(qhg_spinor_field fwd[NS*NC], qhg_spinor_field bwd[NS*NC], qhg_gauge_field gf,
                 int source_coords[ND], qhg_thrp_nn_sink_params thrp_sink)
{  
  qhg_lattice *lat = fwd[0].lat;
  qhg_correlator corr = qhg_correlator_init(SITE_SIZE, lat);
  unsigned long int lvol = lat->lvol;
  unsigned long int **nn = lat->nn;
  for(int i=0; i<ND; i++)
    corr.origin[i] = source_coords[i];
  int Lt = lat->dims[0];
  int tsrc = corr.origin[0];
  int tsnk = (corr.origin[0] + thrp_sink.dt) % Lt;  
  corr.cutoff[0] = tsnk;
  unsigned long int lv3 = lat->lv3;
  int lt = lat->ldims[0];
  int t0 = lat->ldims[0]*lat->comms->proc_coords[0];  

  for(int i=0; i<NS*NC; i++) {
    qhg_xchange_spinor(bwd[i]);
    qhg_xchange_spinor(fwd[i]);
  }
   
  qhg_correlator corr_tmp[4];
  for(int j=0; j<4; j++) {
    corr_tmp[j] = qhg_correlator_init(NCHAN, lat);
  }
  for(int mu=0; mu<ND; mu++) {
    for(int nu=0; nu<ND; nu++) {

      if(lat->comms->proc_id == 0)
        printf("Creating gamma-less local double derivative for mu=%d and nu=%d\n",mu,nu);
#ifdef QHG_OMP
#pragma omp parallel for
#endif
      for(unsigned long int v=0; v<lvol; v++) {

        /* Compute gamma-less local double derivative */
        _Complex double A[4][NS*NC][NS*NC];
        _Complex double T1[NS*NC][NS*NC];
        _Complex double T2[NS*NC][NS*NC];
      
        unsigned long int vp_mu = nn[mu][v];
        unsigned long int vm_mu = nn[mu+ND][v];
        unsigned long int vp_nu = nn[nu][v];
        unsigned long int vm_nu = nn[nu+ND][v];
      
        _Complex double U0_mu[NC][NC];
        _Complex double Um_mu[NC][NC];
        _Complex double U0_nu[NC][NC];
        _Complex double Um_nu[NC][NC];
        su3_load(U0_mu, gf, v, mu);
        su3_load(Um_mu, gf, vm_mu, mu);
        su3_load(U0_nu, gf, v, nu);
        su3_load(Um_nu, gf, vm_nu, nu);
      
        _Complex double FP_mu[NS*NC][NS*NC];
        _Complex double FM_mu[NS*NC][NS*NC];
        _Complex double BP_nu[NS*NC][NS*NC];
        _Complex double BM_nu[NS*NC][NS*NC];
        prop_load(FP_mu, fwd, vp_mu);
        prop_load(FM_mu, fwd, vm_mu);
        prop_load(BP_nu, bwd, vp_nu);
        prop_load(BM_nu, bwd, vm_nu);      
        
        if(t0 + (v/lv3) < tsrc) {
          prop_bc(fwd[0].bc, FP_mu);
          prop_bc(fwd[0].bc, FM_mu);
        }
        
        if(t0 + (v/lv3) > tsnk) {
          prop_bc(bwd[0].bc, BP_nu);
          prop_bc(bwd[0].bc, BM_nu);
        }
        
        if(mu == 0) {
          if(t0 + (v/lv3) == Lt-1 || t0 + (v/lv3) == tsrc-1) {
            /* This catches the case when v is on the right side of the lattice,
               but vp is over the edge or the source */
            prop_bc(fwd[0].bc, FP_mu);
          }	
          if(t0 + (v/lv3) == 0 || t0 + (v/lv3) == tsrc) {
            /* This catches the case when v is on the left side of the lattice,
               but vm is over the edge or the source */
            prop_bc(fwd[0].bc, FM_mu);
          }
        }
        if(nu == 0) {
          if(t0 + (v/lv3) == Lt-1 || t0 + (v/lv3) == tsnk-1) {
            /* This catches the case when v is on the right side of the lattice,
               but vp is over the edge or the sink */
            prop_bc(bwd[0].bc, BP_nu);
          }
          if(t0 + (v/lv3) == 0 || t0 + (v/lv3) == tsnk) {
            /* This catches the case when v is on the left side of the lattice,
               but vm is over the edge or the sink */
            prop_bc(bwd[0].bc, BM_nu);
          }
        }

        // +mu + nu
        prop_mul_su3_U_G(T1, FP_mu, U0_mu);
        prop_mul_su3_D_G(T2, T1, U0_nu);
        prop_mul_gg(A[0], T2, BP_nu);
        
        // -mu - nu
        prop_mul_su3_D_G(T1, FM_mu, Um_mu);
        prop_mul_su3_U_G(T2, T1, Um_nu);
        prop_mul_gg(A[1], T2, BM_nu);
        
        // +mu - nu
        prop_mul_su3_U_G(T1, FP_mu, U0_mu);
        prop_mul_su3_U_G(T2, T1, Um_nu);
        prop_mul_gg(A[2], T2, BM_nu);   
        
        // -mu + nu
        prop_mul_su3_D_G(T1, FM_mu, Um_mu);
        prop_mul_su3_D_G(T2, T1, U0_nu);
        prop_mul_gg(A[3], T2, BP_nu);      

        for(int j=0; j<4; j++) {
          for(int i=0; i<NCHAN; i++) {
            switch(i) {
              /* Scalar */
            case one:
              prop_1_G(T1, A[j]);
              break;
              /* Pseudo-scalar */
            case g5:
              prop_g5_G(T1, A[j]);
              break;
              /* Vector */
            case g0:
              prop_g0_G(T1, A[j]);
              break;
            case gx:
              prop_gx_G(T1, A[j]);
              break;
            case gy:
              prop_gy_G(T1, A[j]);
              break;
            case gz:
              prop_gz_G(T1, A[j]);
              break;
              /* Axial */
            case g5g0:
              prop_g5g0_G(T1, A[j]);
              break;
            case g5gx:
              prop_g5gx_G(T1, A[j]);
              break;
            case g5gy:
              prop_g5gy_G(T1, A[j]);
              break;
            case g5gz:
              prop_g5gz_G(T1, A[j]);
              break;
              /* Tensor */
            case g5si0x:
              prop_g5si0x_G(T1, A[j]);
              break;
            case g5si0y:
              prop_g5si0y_G(T1, A[j]);
              break;
            case g5si0z:
              prop_g5si0z_G(T1, A[j]);
              break;
            case g5sixy:
              prop_g5sixy_G(T1, A[j]);
              break;
            case g5sixz:
              prop_g5sixz_G(T1, A[j]);
              break;
            case g5siyz:
              prop_g5siyz_G(T1, A[j]);
              break;
            }
            // we give minus sign to case 0 and 1
            corr_tmp[j].C[CHAN_IDX(v,i)] = ((j<2)?-1.:+1.)*prop_trace(T1)/16.;
            // adding to correlator
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[j].C[CHAN_IDX(v,i)];
          }//i
        }//j
      }//v

      if(lat->comms->proc_id == 0)
        printf("Created gamma-less local double derivative for mu=%d and nu=%d\n",mu,nu);
      
      // shifting and adding to correlator
      int shift[ND] = {0};

      if(mu!=nu) {
        // +mu
        shift[mu]=-1;
        shift[nu]=0;
        qhg_correlator_shift(corr_tmp[0], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[0].C[CHAN_IDX(v,i)];
          }
        }
        
        // +mu +nu
        shift[mu]=0;
        shift[nu]=-1;
        qhg_correlator_shift(corr_tmp[0], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[0].C[CHAN_IDX(v,i)];
          }
        }
        
        // +nu
        shift[mu]=+1;
        shift[nu]=0;
        qhg_correlator_shift(corr_tmp[0], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[0].C[CHAN_IDX(v,i)];
          }
        }
        
        // -mu
        shift[mu]=+1;
        shift[nu]=0;
        qhg_correlator_shift(corr_tmp[1], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[1].C[CHAN_IDX(v,i)];
          }
        }
        
        // -mu -nu
        shift[mu]=0;
        shift[nu]=+1;
        qhg_correlator_shift(corr_tmp[1], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[1].C[CHAN_IDX(v,i)];
          }
        }
        
        // -nu
        shift[mu]=-1;
        shift[nu]=0;
        qhg_correlator_shift(corr_tmp[1], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[1].C[CHAN_IDX(v,i)];
          }
        }      
        
        // +mu
        shift[mu]=-1;
        shift[nu]=0;
        qhg_correlator_shift(corr_tmp[2], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[2].C[CHAN_IDX(v,i)];
          }
        }
        
        // +mu -nu
        shift[mu]=0;
        shift[nu]=+1;
        qhg_correlator_shift(corr_tmp[2], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[2].C[CHAN_IDX(v,i)];
          }
        }
        
        // -nu
        shift[mu]=+1;
        shift[nu]=0;
        qhg_correlator_shift(corr_tmp[2], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[2].C[CHAN_IDX(v,i)];
          }
        }      
        
        // -mu
        shift[mu]=+1;
        shift[nu]=0;
        qhg_correlator_shift(corr_tmp[3], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[3].C[CHAN_IDX(v,i)];
          }
        }
        
        // -mu +nu
        shift[mu]=0;
        shift[nu]=-1;
        qhg_correlator_shift(corr_tmp[3], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[3].C[CHAN_IDX(v,i)];
          }
        }
        
        // +nu
        shift[mu]=-1;
        shift[nu]=0;
        qhg_correlator_shift(corr_tmp[3], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[3].C[CHAN_IDX(v,i)];
          }
        }
      }
      else {
        // +mu
        shift[mu]=-1;
        qhg_correlator_shift(corr_tmp[0], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += 2*corr_tmp[0].C[CHAN_IDX(v,i)];
          }
        }
        
        // +2*mu
        shift[mu]=-1;
        qhg_correlator_shift(corr_tmp[0], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[0].C[CHAN_IDX(v,i)];
          }
        }
        
        // -mu
        shift[mu]=+1;
        qhg_correlator_shift(corr_tmp[1], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += 2*corr_tmp[1].C[CHAN_IDX(v,i)];
          }
        }
        
        // -2*mu
        shift[mu]=+1;
        qhg_correlator_shift(corr_tmp[1], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[1].C[CHAN_IDX(v,i)];
          }
        }

        // +mu
        shift[mu]=-1;
        qhg_correlator_shift(corr_tmp[2], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[2].C[CHAN_IDX(v,i)];
          }
        }
        
        // 0
        shift[mu]=+1;
        qhg_correlator_shift(corr_tmp[2], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[2].C[CHAN_IDX(v,i)];
          }
        }
        
        // -mu
        shift[mu]=+1;
        qhg_correlator_shift(corr_tmp[2], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[2].C[CHAN_IDX(v,i)];
          }
        }      
        
        // -mu
        shift[mu]=+1;
        qhg_correlator_shift(corr_tmp[3], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[3].C[CHAN_IDX(v,i)];
          }
        }
        
        // 0
        shift[mu]=-1;
        qhg_correlator_shift(corr_tmp[3], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[3].C[CHAN_IDX(v,i)];
          }
        }
        
        // +mu
        shift[mu]=-1;
        qhg_correlator_shift(corr_tmp[3], shift);
        for(unsigned long int v=0; v<lvol; v++) {
          for(int i=0; i<NCHAN; i++) {
            corr.C[TOT_IDX(v,mu,nu,i)] += corr_tmp[3].C[CHAN_IDX(v,i)];
          }
        }

      }
      if(lat->comms->proc_id == 0)
        printf("correlator completed for mu=%d and nu=%d\n",mu,nu);
    }//nu
  }//mu

  for(int j=0; j<4; j++) {
    qhg_correlator_finalize(corr_tmp[j]);
  }
    
  corr.mom_list = NULL;
  return corr;
}

qhg_correlator
qhg_nn_thrp_der2_test(qhg_spinor_field fwd[NS*NC], qhg_spinor_field bwd[NS*NC], qhg_gauge_field gf,
                      int source_coords[ND], qhg_thrp_nn_sink_params thrp_sink)
{  
  qhg_lattice *lat = fwd[0].lat;
  qhg_correlator corr = qhg_correlator_init(SITE_SIZE, lat);
  unsigned long int lvol = lat->lvol;
  unsigned long int **nn = lat->nn;
  for(int i=0; i<ND; i++)
    corr.origin[i] = source_coords[i];
  int Lt = lat->dims[0];
  int tsrc = corr.origin[0];
  int tsnk = (corr.origin[0] + thrp_sink.dt) % Lt;  
  corr.cutoff[0] = tsnk;
  unsigned long int lv3 = lat->lv3;
  int lt = lat->ldims[0];
  int t0 = lat->ldims[0]*lat->comms->proc_coords[0];  

  for(int i=0; i<NS*NC; i++) {
    qhg_xchange_spinor(bwd[i]);
    qhg_xchange_spinor(fwd[i]);
  }

  for(int mu=0; mu<ND; mu++) {
    for(int nu=0; nu<ND; nu++) {
      if( mu == nu) continue;
      if(lat->comms->proc_id == 0)
        printf("Creating double derivative for mu=%d and nu=%d\n",mu,nu);

      for(unsigned long int v=0; v<lvol; v++) {

        /* Compute gamma-less double derivative */
        _Complex double A[NS*NC][NS*NC]; //accumulation
        _Complex double T1[NS*NC][NS*NC]; //temporal
        _Complex double T2[NS*NC][NS*NC]; //temporal
        for( int i=0; i<NS*NC; i++)
          for( int j=0; j<NS*NC; j++) {
            A[i][j]=0;
          }
      
        int dir[2] = {mu,nu};
        unsigned long int v_m0p[3][3];
        v_m0p[1][0] = nn[dir[1]+ND][v]; //minus nu
        v_m0p[1][1] = v; //0
        v_m0p[1][2] = nn[dir[1]][v]; //plus nu
        for( int i=0; i<3; i++) {
          v_m0p[0][i]= nn[dir[0]+ND][v_m0p[1][i]]; 
          v_m0p[2][i]= nn[dir[0]][v_m0p[1][i]]; 
        }
      
        _Complex double U[2][3][2][NC][NC];
        for(int i=0; i<3; i++) { // U[1] links in dir 0
          su3_load_D(U[0][i][0], gf, v_m0p[0][i], dir[0]); //backward
          su3_load(U[0][i][1], gf, v_m0p[1][i], dir[0]); // forward
        }
        for(int i=0; i<3; i++) { // U[1] links in dir 1
          su3_load_D(U[1][i][0], gf, v_m0p[i][0], dir[1]); //backward
          su3_load(U[1][i][1], gf, v_m0p[i][1], dir[1]); // forward
        }
      
        _Complex double F[3][3][NS*NC][NS*NC];
        _Complex double B[3][3][NS*NC][NS*NC];
        for( int i=0; i<3; i++)
          for( int j=0; j<3; j++) {
            prop_load(F[i][j], fwd, v_m0p[i][j]);
            prop_load(B[i][j], bwd, v_m0p[i][j]);
          }
      
        if(t0 + (v/lv3) < tsrc)
          for( int i=0; i<3; i++)
            for( int j=0; j<3; j++) {
              prop_bc(fwd[0].bc, F[i][j]);
            }
        
        if(t0 + (v/lv3) > tsnk)
          for( int i=0; i<3; i++)
            for( int j=0; j<3; j++) {
              prop_bc(bwd[0].bc, B[i][j]);
            }

        if(dir[0] == 0) {
          if(t0 + (v/lv3) == Lt-1 || t0 + (v/lv3) == tsrc-1) {
            /* This catches the case when v is on the right side of the lattice,
               but vp is over the edge or the source */
            for( int i=0; i<3; i++) {
              prop_bc(fwd[0].bc, F[2][i]); 
            }
          }	
          if(t0 + (v/lv3) == Lt-1 || t0 + (v/lv3) == tsnk) {
            /* This catches the case when v is on the right side of the lattice,
               but vp is over the edge or the sink */
            for( int i=0; i<3; i++) {
              prop_bc(bwd[0].bc, B[2][i]); 
            }
          }	
          if(t0 + (v/lv3) == 0 || t0 + (v/lv3) == tsrc) {
            /* This catches the case when v is on the right side of the lattice,
               but vm is over the edge or the source */
            for( int i=0; i<3; i++) {
              prop_bc(fwd[0].bc, F[0][i]); 
            }
          }	
          if(t0 + (v/lv3) == 0 || t0 + (v/lv3) == tsnk + 1) {
            /* This catches the case when v is on the right side of the lattice,
               but vm is over the edge or the sink */
            for( int i=0; i<3; i++) {
              prop_bc(fwd[0].bc, B[0][i]); 
            }
          }	
        }
        if(dir[1] == 0) {
          if(t0 + (v/lv3) == Lt-1 || t0 + (v/lv3) == tsrc-1) {
            /* This catches the case when v is on the right side of the lattice,
               but vp is over the edge or the source */
            for( int i=0; i<3; i++) {
              prop_bc(fwd[0].bc, F[i][2]); 
            }
          }	
          if(t0 + (v/lv3) == Lt-1 || t0 + (v/lv3) == tsnk) {
            /* This catches the case when v is on the right side of the lattice,
               but vp is over the edge or the sink */
            for( int i=0; i<3; i++) {
              prop_bc(bwd[0].bc, B[i][2]); 
            }
          }	
          if(t0 + (v/lv3) == 0 || t0 + (v/lv3) == tsrc) {
            /* This catches the case when v is on the right side of the lattice,
               but vm is over the edge or the source */
            for( int i=0; i<3; i++) {
              prop_bc(fwd[0].bc, F[i][0]); 
            }
          }	
          if(t0 + (v/lv3) == 0 || t0 + (v/lv3) == tsnk + 1) {
            /* This catches the case when v is on the right side of the lattice,
               but vm is over the edge or the sink */
            for( int i=0; i<3; i++) {
              prop_bc(fwd[0].bc, B[i][0]); 
            }
          }	
        }

        // DD(-> ->) goes from 4 corner points to the center
        for( int i=0; i<3; i+=2) 
          for( int j=0; j<3; j+=2) {
            prop_mul_su3_U_G(T1, F[i][j], U[0][j][i==0?0:1]);
            prop_mul_su3_U_G(T2, T1, U[1][1][j==0?0:1]);
            prop_mul_gg(T1, T2, B[1][1]);
            if(i==j)
              prop_peq_g(A,T1);
            else
              prop_meq_g(A,T1);
          }

        // DD(<- <-) goes the center to the 4 corner points
        for( int i=0; i<3; i+=2) 
          for( int j=0; j<3; j+=2) {
            prop_mul_su3_D_G(T1, F[1][1], U[0][1][i==0?0:1]);
            prop_mul_su3_D_G(T2, T1, U[1][i][j==0?0:1]);
            prop_mul_gg(T1, T2, B[i][j]);
            if(i==j)
              prop_peq_g(A,T1);
            else
              prop_meq_g(A,T1);
          }

        // DD(<- ->) goes from up/dn to lf/rg through the center
        for( int i=0; i<3; i+=2) 
          for( int j=0; j<3; j+=2) {
            prop_mul_su3_U_G(T1, F[i][1], U[0][1][i==0?0:1]);
            prop_mul_su3_D_G(T2, T1, U[1][1][j==0?0:1]);
            prop_mul_gg(T1, T2, B[1][j]);
            if(i==j)
              prop_meq_g(A,T1);
            else
              prop_peq_g(A,T1);
          }

        // DD(-> <-) goes from lf/rg to up/dn around
        for( int i=0; i<3; i+=2) 
          for( int j=0; j<3; j+=2) {
            prop_mul_su3_D_G(T1, F[1][j], U[0][j][i==0?0:1]);
            prop_mul_su3_U_G(T2, T1, U[1][i][j==0?0:1]);
            prop_mul_gg(T1, T2, B[i][1]);
            if(i==j)
              prop_meq_g(A,T1);
            else
              prop_peq_g(A,T1);
          }

        for(int i=0; i<NCHAN; i++) {
          switch(i) {
            /* Scalar */
          case one:
            prop_1_G(T1, A);
            break;
            /* Pseudo-scalar */
          case g5:
            prop_g5_G(T1, A);
            break;
            /* Vector */
          case g0:
            prop_g0_G(T1, A);
            break;
          case gx:
            prop_gx_G(T1, A);
            break;
          case gy:
            prop_gy_G(T1, A);
            break;
          case gz:
            prop_gz_G(T1, A);
            break;
            /* Axial */
          case g5g0:
            prop_g5g0_G(T1, A);
            break;
          case g5gx:
            prop_g5gx_G(T1, A);
            break;
          case g5gy:
            prop_g5gy_G(T1, A);
            break;
          case g5gz:
            prop_g5gz_G(T1, A);
            break;
            /* Tensor */
          case g5si0x:
            prop_g5si0x_G(T1, A);
            break;
          case g5si0y:
            prop_g5si0y_G(T1, A);
            break;
          case g5si0z:
            prop_g5si0z_G(T1, A);
            break;
          case g5sixy:
            prop_g5sixy_G(T1, A);
            break;
            case g5sixz:
              prop_g5sixz_G(T1, A);
              break;
          case g5siyz:
            prop_g5siyz_G(T1, A);
            break;
          }
          // saving into correlator
          corr.C[TOT_IDX(v,mu,nu,i)] = prop_trace(T1)/16.;
        }//i
      }//v

      if(lat->comms->proc_id == 0)
        printf("Created double derivative for mu=%d and nu=%d\n",mu,nu);
    }//nu
  }//mu

  corr.mom_list = NULL;
  return corr;
}
