#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include <complex.h>
#include <qhg_idx.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_correlator.h>
#include <qhg_correlator_shift.h>
#include <qhg_spinor_field.h>
#include <qhg_fast_spinor_field.h>
#include <qhg_fast_gauge_field.h>
#include <qhg_xchange_spinor.h>
#include <qhg_xchange_gauge.h>
#include <qhg_prop_gammas.h>
#include <qhg_prop_ops.h>
#include <qhg_su3_ops.h>
#include <qhg_io_utils.h>
#include <qhg_nn_thrp_der.h>
#include <math.h>

unsigned int fact(unsigned int n) {
  if(n > 1)
    return n*fact(n-1);
  else
    return 1;
}

void swap(int *x, int *y) { 
  int temp; 
  temp = *x; 
  *x = *y; 
  *y = temp; 
} 

void permute(int *dirs, int l, int r, int **save) { 
  if (l == r-1) { 
    for(int i=0; i<r; i++) {
      (*save)[i] = dirs[i];
    }
    *save += r; 
  }
  else { 
    for (int i = l; i < r; i++) { 
      if(i==l || (dirs[l])!=(dirs[i])) {
        swap((dirs+l), (dirs+i)); 
        permute(dirs, l+1, r, save); 
        swap((dirs+l), (dirs+i)); 
      }
    } 
  } 
} 

const int gamma_idx[16][4] = {
  { 0, 5,10,15}, //1
  { 0, 5,10,15}, //g5
  { 2, 7, 8,13}, //g0
  { 3, 6, 9,12}, //gx
  { 3, 6, 9,12}, //gy
  { 2, 7, 8,13}, //gz
  { 2, 7, 8,13}, //g5g0
  { 3, 6, 9,12}, //g5gx
  { 3, 6, 9,12}, //g5gy
  { 2, 7, 8,13}, //g5gz
  { 1, 4,11,14}, //g5si0x
  { 1, 4,11,14}, //g5si0y
  { 0, 5,10,15}, //g5si0z
  { 0, 5,10,15}, //g5sixy
  { 1, 4,11,14}, //g5sixz
  { 1, 4,11,14}, //g5siyz
};

const _Complex double gamma_value[16][4] = {
  {+1,+1,+1,+1}, //1
  {+1,+1,-1,-1}, //g5
  {-1,-1,-1,-1}, //g0
  {+I,+I,-I,-I}, //gx
  {-1,+1,+1,-1}, //gy
  {+I,-I,-I,+I}, //gz
  {+1,+1,-1,-1}, //g5g0
  {-I,-I,-I,-I}, //g5gx
  {+1,-1,+1,-1}, //g5gy
  {-I,+I,-I,+I}, //g5gz
  {+I,+I,+I,+I}, //g5si0x
  {-1,+1,-1,+1}, //g5si0y
  {+I,-I,+I,-I}, //g5si0z
  {-I,+I,+I,-I}, //g5sixy
  {-1,+1,+1,-1}, //g5sixz
  {-I,-I,+I,+I}, //g5siyz
};

qhg_der_correlator
qhg_nn_thrp_der(qhg_spinor_field fwd[NS*NC], qhg_spinor_field bwd[NS*NC], qhg_gauge_field gf,
                int source_coords[ND], qhg_thrp_nn_sink_params thrp_sink, int der_order, bool (*to_skip)(int*), bool rev_trick, bool perm_trick)
{
  qhg_lattice *lat = gf.lat;
  qhg_comms *comms = lat->comms;
  int Lt = lat->dims[0];
  qhg_fast_gauge_field gf_1, gf_2, gf_3;
  qhg_fast_spinor_field sfwd, sbwd, shifted_sfwd, shifted_sbwd, stmp;
  qhg_der_correlator corr;
  qhg_correlator corr_tmp;
  qhg_mom_list *mom_list;
  int nmoms;
  int (*moms)[3];
  int site_size = NS*NS;
  size_t vol_size;
  corr.lat = lat;
  corr.der_order = der_order;
  corr.site_size = site_size;
  corr.dt = thrp_sink.dt;
  corr.proj = thrp_sink.proj;
  corr.ncorr = 1<<(2*der_order);
  corr.origin = malloc(4*sizeof(int));
  for(int i=0; i < 4; i++)
    corr.origin[i] = source_coords[i];
  if( thrp_sink.corr_space == MOM_SPACE ) {
    corr.mom_list = thrp_sink.mom_list;
    nmoms = corr.mom_list->n_mom_vecs;
    moms = corr.mom_list->mom_vecs;
    vol_size = lat->dims[0]*nmoms;
  } else {
    vol_size = lat->lvol;
    corr.mom_list = NULL;
  }
  corr.C = qhg_alloc(corr.ncorr*sizeof(_Complex double *));
  for(int i=0; i<corr.ncorr; i++) { 
    corr.C[i]=NULL;
  }

  // allocation
  sbwd = qhg_fast_spinor_field_init(bwd[0].lat, bwd[0].bc);  
  sfwd = qhg_fast_spinor_field_init(fwd[0].lat, fwd[0].bc);  
  shifted_sbwd = qhg_fast_spinor_field_init(bwd[0].lat, bwd[0].bc);  
  shifted_sfwd = qhg_fast_spinor_field_init(fwd[0].lat, fwd[0].bc);  
  stmp = qhg_fast_spinor_field_init(fwd[0].lat, fwd[0].bc);  
  gf_1 = qhg_fast_gauge_field_init(lat);
  if( der_order > 1 ) {
    gf_2 = qhg_fast_gauge_field_init(lat);
    gf_3 = qhg_fast_gauge_field_init(lat);
  }

  // loading bwd and removing boundary conditions
  qhg_fast_spinor_field_import(sbwd, bwd);
  qhg_fast_spinor_field_remove_bc_higher(sbwd, (source_coords[0] + thrp_sink.dt ) % Lt);
  qhg_fast_spinor_field_import(sfwd, fwd);
  qhg_fast_spinor_field_remove_bc_lower(sfwd, source_coords[0]);

  // loop over all dirs (4) and signs (2) combinations: (4*2)^der_order
  for(int i=0; i<(1<<(3*der_order)); i++) { 

    // extracting information from i: directions, signs, shift, corr_id
    int tmp_i = i;
    int dirs[der_order];
    int end[ND] = {0,0,0,0};
    int count_dirs[2*ND] = {0,0,0,0,0,0,0,0}; 
    double factor = 1./pow(2,2*der_order);
    for(int j=0; j< der_order; j++) {
      dirs[j] = tmp_i % 8;
      tmp_i /= 8;
      count_dirs[dirs[j]]++;
      end[dirs[j]%ND] += (1 - 2*(dirs[j]/ND));
      factor *= (1. - 2.*(dirs[j]/ND)); 
    }
    if(to_skip != NULL && to_skip(count_dirs)) continue;
    
    int n_perms = 1;
    if(perm_trick) {
      // Counting the number of permutations without repeated elements
      n_perms = fact(der_order);
      for(int j=0; j< 2*ND; j++) {
        n_perms /= fact(count_dirs[j]);
      }
    }
    int perm_dirs[n_perms][der_order];
    if(perm_trick) {
      // List of all the permutations of the directions
      int *save = perm_dirs[0];
      permute(dirs, 0, der_order, &save);
    } else {
      for(int j=0; j<der_order; j++) {
        perm_dirs[0][j] = dirs[j];
      }
    }

    int max_corr = (rev_trick ? 2 : 1) * n_perms;
    int id[max_corr];
    int corr_id[max_corr];
    bool is_dag[max_corr];
    int n_corr=0;
    bool do_rev=false;
    bool skip = false;
    for(int j=0; j<n_perms; j++) {
      id[n_corr] = 0;
      corr_id[n_corr] = 0;
      is_dag[n_corr] = false;
      for(int k=0; k<der_order; k++) {
        id[n_corr] = id[n_corr]*2*ND + perm_dirs[j][k];
        corr_id[n_corr] = corr_id[n_corr]*ND + perm_dirs[j][k]%ND;        
      }
      if(id[n_corr] < id[0]) { skip=true; break; } // we have already done it
      n_corr++;
      if(rev_trick) {
        id[n_corr] = 0;
        corr_id[n_corr] = 0;
        is_dag[n_corr] = true;
        for(int k=0; k<der_order; k++) {
          id[n_corr] = id[n_corr]*2*ND + (perm_dirs[j][der_order-1-k] + ND)%(2*ND);
          corr_id[n_corr] = corr_id[n_corr]*ND + perm_dirs[j][der_order-1-k]%ND;        
        }
        if(id[n_corr] < id[0]) { // we have already done it 
          skip=true; break; 
        } else if(id[n_corr] != id[n_corr-1]){ // we don't want to do the same id twice
          n_corr++; do_rev=true; 
        } 
      }
    }
    if(skip) {
      continue;
    }

    for(int j=0; j<n_corr; j++) {
      if(corr.C[corr_id[j]] == NULL) {
        corr.C[corr_id[j]] = qhg_alloc(vol_size*site_size*sizeof(_Complex double));
        memset(corr.C[corr_id[j]], '\0', vol_size*site_size*sizeof(_Complex double));
      }
    }
    
    // shifting the propagator
    bool shifted = false;
    for(int dir=0; dir < ND; dir++) {
      for(int j=0; j<abs(end[dir]); j++) {
        qhg_fast_spinor_field_shift(stmp, shifted ? shifted_sfwd : sfwd, dir + ((end[dir]>0) ? 0:ND));
        qhg_fast_spinor_field_swap(&shifted_sfwd, &stmp);
        if(do_rev) {
          qhg_fast_spinor_field_shift(stmp, shifted ? shifted_sbwd : sbwd, dir + ((end[dir]>0) ? 0:ND));
          qhg_fast_spinor_field_swap(&shifted_sbwd, &stmp);
        }
        shifted = true;
      }
    }
    // in case we didn't shift at all we still need to copy
    if(shifted == false) {
      qhg_fast_spinor_field_copy(shifted_sfwd, sfwd);
      if(do_rev) qhg_fast_spinor_field_copy(shifted_sbwd, sbwd);
    }
    // contracting the propagators
    qhg_fast_spinor_field_multiply_G_G(stmp, shifted_sfwd, sbwd);
    qhg_fast_spinor_field_swap(&shifted_sfwd, &stmp);
    if(do_rev) {
      qhg_fast_spinor_field_multiply_G_G(stmp, sfwd, shifted_sbwd);
      qhg_fast_spinor_field_swap(&shifted_sbwd, &stmp);
    }

    tmp_i = 0;
    qhg_correlator corr_tmp = qhg_correlator_init(n_corr*NS*NS, lat);
    for(int p=0; p<n_perms; p++) {
      // constructing the path in the reverse order
      for(int j=der_order-1; j >= 0; j--) {
        // skipping trivial multiplications
        if(j>0 && perm_dirs[p][j] == (perm_dirs[p][j-1]+ND)%(2*ND) ) {
          if(j==der_order-1) qhg_fast_gauge_field_id(gf_1); // initializing to id
          j--;
          continue;
        }
        if ( j==der_order-1 ) {
          // loading the first dir from gf and saving in gf_1
          qhg_fast_gauge_field_import_dir(gf_1, gf, perm_dirs[p][j]);
        } else {
          qhg_fast_gauge_field_shift(gf_3, gf_1, perm_dirs[p][j]);
          qhg_fast_gauge_field_import_dir(gf_2, gf, perm_dirs[p][j]);
          qhg_fast_gauge_field_multiply(gf_1, gf_2, gf_3);
        }
      }

      // Closing the path
      qhg_fast_spinor_field_trace_multiply_U_G( gf_1, shifted_sfwd, corr_tmp, tmp_i, factor);
      tmp_i++;
      if(rev_trick && tmp_i < n_corr && is_dag[tmp_i]) {
        qhg_fast_spinor_field_trace_multiply_Udag_G( gf_1, shifted_sbwd, corr_tmp, tmp_i,
                                                     (der_order%2==1) ? -factor:factor);
        tmp_i++;
      }
    }

    // Position space or momentum space
    if(thrp_sink.corr_space == POS_SPACE) {
      /*
       * The same term needs to be summed in all the sites of the kind
       * Order 1: 0, +A
       * Order 2: 0, +A, +A+B, +B
       * Order 3: 0, +A, +A+B, +B, +B+C, +A+B+C, +A+C, +C
       * where A, B, C, etc are the different signed directions.
       * The list is constructed by adding a new direction and reading backwards the previous list.
       */
      int add_to[1<<der_order][4]; //2^der_order components, 4 dirs
      for(int d=0; d<4; d++) {
        add_to[0][d] = 0;
      }
      for(int j=0; j<der_order; j++) {
        for(int k=(1<<j); k>0; k--) {
          for(int d=0; d<4; d++) {
            add_to[(1<<(j+1))-k][d] = add_to[k-1][d];
          }
          add_to[(1<<(j+1))-k][dirs[j]%ND] -= (1-2*(dirs[j]/ND)); 
        } 
      }
      // looping over all the add_to: shifting and adding
      for(int j=0; j<(1<<der_order); j++) {
        int diff[4];
        for(int d=0; d<4; d++) {
          diff[d] = j==0 ? (add_to[0][d]) : (add_to[j][d] - add_to[j-1][d]);
        }
        qhg_correlator_shift(corr_tmp, diff);
        for(size_t v=0; v<vol_size; v++) {
          for(int k=0; k<n_corr; k++) {
            for(int s=0; s<NS*NS; s++) {
              corr.C[corr_id[k]][v*NS*NS+s] += corr_tmp.C[(v*n_corr+k)*NS*NS + s];
            }
          }
        }
      }
    }
    else if( thrp_sink.corr_space == MOM_SPACE ) {
      /*
       * As above:
       * The same term needs to be summed in all the sites of the kind
       * Order 1: 0, +A
       * Order 2: 0, +A, +A+B, +B
       * Order 3: 0, +A, +A+B, +B, +B+C, +A+B+C, +A+C, +C
       * where A, B, C, etc are the different signed directions.
       * The list is constructed by adding a new direction and reading backwards the previous list.
       *
       * Here we split between time and space contributions. 
       * The space contributions are taken into account in the fourier transform.
       * The time contributions are added after by shifting in time.
       */
      int n_time = 0;
      for(int j=0; j<der_order; j++)
        if(dirs[j]%ND==0) n_time++;
      int n_space = der_order-n_time;
      int add_to_time[1<<n_time]; //2^n_time components, 1 dirs
      int add_to_space[1<<n_space][3]; //2^n_space components, 3 dirs
      add_to_time[0] = 0;
      for(int d=0; d<3; d++)
        add_to_space[0][d] = 0;

      int it = 0, is =0;
      for(int j=0; j<der_order; j++) {
        if(dirs[j]%ND == 0) {
          for(int k=(1<<it); k>0; k--) {
            add_to_time[(1<<(it+1))-k] = add_to_time[k-1] + (1-2*(dirs[j]/ND));
          } 
          it++;
        } else {
          for(int k=(1<<is); k>0; k--) {
            for(int d=0; d<3; d++)
              add_to_space[(1<<(is+1))-k][d] = add_to_space[k-1][d];
            add_to_space[(1<<(is+1))-k][dirs[j]%ND-1] += (1-2*(dirs[j]/ND)); 
          }
          is++;
        }
      }

      // computing the shift for going to global coordinates centered in the source
      double shift[3];
      for(int d=0; d<3; d++)
        shift[d] = lat->dims[d+1] + lat->ldims[d+1]*comms->proc_coords[d+1] - source_coords[d+1];

      _Complex double accum[nmoms*lat->ldims[0]*site_size*n_corr];
      // Fourier transform
      for(int m=0; m<nmoms; m++) {
        int *k = moms[m];
        // Global coefficient due to the add_to_space terms
        _Complex double coeff = 0;
        for(int j=0; j<(1<<(n_space)); j++) {
          double ph = 0;
          for(int d=0; d<3; d++)
            ph += (k[d]*add_to_space[j][d])/lat->dims[d+1];
          ph *= 2*M_PI;
          coeff += cos(ph) - _Complex_I * sin(ph);
        }
        // Sum over all space
        int id[4];
        for(id[1]=0; id[1]<lat->ldims[1]; id[1]++)
          for(id[2]=0; id[2]<lat->ldims[2]; id[2]++)
            for(id[3]=0; id[3]<lat->ldims[3]; id[3]++) {
              double ph=0;
              for(int d=0; d<3; d++)
                ph += k[d]*(id[d+1]+shift[d])/lat->dims[d+1];
              ph *= 2*M_PI;
              _Complex double loc_coeff = cos(ph) - _Complex_I * sin(ph);
              loc_coeff *= coeff;
              for(id[0]=0; id[0]<lat->ldims[0]; id[0]++) {
                unsigned long int v = IDX(id, lat->ldims);
                for(int s=0; s<site_size*n_corr; s++) {
                  accum[(id[0]*nmoms+m)*site_size*n_corr + s] += loc_coeff*corr_tmp.C[v*site_size*n_corr + s];
                }
              }
            }
      }

      MPI_Comm space_comm, time_comm;
      int space_rank, time_rank;
      MPI_Comm_split(comms->comm, comms->proc_coords[0], comms->proc_id, &space_comm);
      MPI_Comm_rank( space_comm, &space_rank );
      int pid = (comms->proc_coords[1]*comms->proc_dims[2] + comms->proc_coords[2])*comms->proc_dims[3] + comms->proc_coords[3];
      MPI_Comm_split(comms->comm, pid, comms->proc_coords[0], &time_comm);
      MPI_Comm_rank( time_comm, &time_rank );

      qhg_correlator corr_ft = qhg_transformed_correlator_init(site_size*n_corr, lat, thrp_sink.mom_list);
      MPI_Reduce(accum, corr_ft.C, 2*nmoms*lat->ldims[0]*site_size*n_corr, MPI_DOUBLE, MPI_SUM, 0, space_comm);

      if(pid==0) {
        for(int j=0; j<(1 << n_time); j++) {
          int shift = j==0 ? add_to_time[0] : (add_to_time[j] - add_to_time[j-1]);
          qhg_correlator_shift_ft(corr_ft, shift, time_comm);
          for(size_t v=0; v<vol_size; v++) {
            for(int k=0; k<n_corr; k++) {
              for(int s=0; s<NS*NS; s++) {
                corr.C[corr_id[k]][v*NS*NS+s] += corr_tmp.C[(v*n_corr+k)*NS*NS + s];
              }
            }
          }
        }
      }
      qhg_correlator_finalize(corr_ft);
    }
    else {
      printf("ERROR wrong corr_space\n");
    }
    qhg_correlator_finalize(corr_tmp); 
  }
  for(int k=0; k<corr.ncorr; k++)
    for(size_t v=0; v<vol_size; v++) {
      _Complex double C[NS*NS];
      for(int s=0; s<NS*NS; s++) {
        C[s] = corr.C[k][v*NS*NS+s];
        corr.C[k][v*NS*NS+s] = 0;
      }
      for(int s=0; s<NS*NS; s++)
        for(int j=0; j<4; j++)
          corr.C[k][v*NS*NS+s] += gamma_value[s][j]*C[gamma_idx[s][j]];
    }
  // finalize
  qhg_fast_spinor_field_finalize(sbwd);  
  qhg_fast_spinor_field_finalize(sfwd);  
  qhg_fast_spinor_field_finalize(shifted_sbwd);  
  qhg_fast_spinor_field_finalize(shifted_sfwd);  
  qhg_fast_spinor_field_finalize(stmp);  
  qhg_fast_gauge_field_finalize(gf_1);
  if( der_order > 1 ) {
    qhg_fast_gauge_field_finalize(gf_2);
    qhg_fast_gauge_field_finalize(gf_3);
  }

  return corr;
}
