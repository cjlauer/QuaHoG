#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <unistd.h>

#include <qhg_idx.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_io_utils.h>

int QHG_EPS[NEPS][NC+1] = {{0,1,2,+1},
                           {0,2,1,-1},
                           {1,0,2,-1},
                           {1,2,0,+1},
                           {2,0,1,+1},
                           {2,1,0,-1}};

inline void
qhg_fast_contract_f111(qhg_baryons_open_correlator corr, qhg_fast_spinor_field sp_1) 
{

  qhg_lattice *lat = sp_1.lat; 
  int lt = lat->ldims[0];
  unsigned long int lv3 = lat->lv3;

#ifdef QHG_OMP
#pragma omp parallel for
#endif
  for(int t=0; t<lt; t++) {
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++)
        for(int s2=0; s2<NS; s2++)
          for(int s3=0; s3<NS; s3++)
            for(int s4=0; s4<NS; s4++)
              for(int s5=0; s5<NS; s5++)
                for(int eps0=0; eps0<NEPS; eps0++) {
                  int c4 = QHG_EPS[eps0][0];
                  int c5 = QHG_EPS[eps0][1];
                  int c0 = QHG_EPS[eps0][2];
                  for(int eps1=0; eps1<NEPS; eps1++) {
                    int c2 = QHG_EPS[eps1][0];
                    int c3 = QHG_EPS[eps1][1];
                    int c1 = QHG_EPS[eps1][2];
                    double sign = QHG_EPS[eps0][3]*QHG_EPS[eps1][3];

                    /*
                     * NOTE: the loops have been packed in this way in order to induce self-vectorization
                     */
                    // [10] ( [25][34] - [24][35] )
                    afloat par10_re[lv3];
                    afloat par10_im[lv3];
                    // [14] ( [20][35] - [25][30])
                    afloat par14_re[lv3];
                    afloat par14_im[lv3];
                    // [15] ( [24][30] - [20][34] )
                    afloat par15_re[lv3];
                    afloat par15_im[lv3];

                    // list of restrict pointers
                    afloat * restrict f10_re = sp_1.field[t][s1][s0][c1][c0][0];
                    afloat * restrict f10_im = sp_1.field[t][s1][s0][c1][c0][1];
                    afloat * restrict f14_re = sp_1.field[t][s1][s4][c1][c4][0];
                    afloat * restrict f14_im = sp_1.field[t][s1][s4][c1][c4][1];
                    afloat * restrict f15_re = sp_1.field[t][s1][s5][c1][c5][0];
                    afloat * restrict f15_im = sp_1.field[t][s1][s5][c1][c5][1];

                    afloat * restrict f20_re = sp_1.field[t][s2][s0][c2][c0][0];
                    afloat * restrict f20_im = sp_1.field[t][s2][s0][c2][c0][1];
                    afloat * restrict f24_re = sp_1.field[t][s2][s4][c2][c4][0];
                    afloat * restrict f24_im = sp_1.field[t][s2][s4][c2][c4][1];
                    afloat * restrict f25_re = sp_1.field[t][s2][s5][c2][c5][0];
                    afloat * restrict f25_im = sp_1.field[t][s2][s5][c2][c5][1];

                    afloat * restrict f30_re = sp_1.field[t][s3][s0][c3][c0][0];
                    afloat * restrict f30_im = sp_1.field[t][s3][s0][c3][c0][1];
                    afloat * restrict f34_re = sp_1.field[t][s3][s4][c3][c4][0];
                    afloat * restrict f34_im = sp_1.field[t][s3][s4][c3][c4][1];
                    afloat * restrict f35_re = sp_1.field[t][s3][s5][c3][c5][0];
                    afloat * restrict f35_im = sp_1.field[t][s3][s5][c3][c5][1];

                    for(unsigned long int v=0; v<lv3; v++) {
                      // +[25][34]
                      par10_re[v]  = + ( f25_re[v] * f34_re[v] - f25_im[v] * f34_im[v] );
                      par10_im[v]  = + ( f25_re[v] * f34_im[v] + f25_im[v] * f34_re[v] );
                      // -[25][30]
                      par14_re[v]  = - ( f25_re[v] * f30_re[v] - f25_im[v] * f30_im[v] );
                      par14_im[v]  = - ( f25_re[v] * f30_im[v] + f25_im[v] * f30_re[v] );
                      // +[24][30]
                      par15_re[v]  = + ( f24_re[v] * f30_re[v] - f24_im[v] * f30_im[v] );
                      par15_im[v]  = + ( f24_re[v] * f30_im[v] + f24_im[v] * f30_re[v] );
                      // -[24][35]
                      par10_re[v] += - ( f24_re[v] * f35_re[v] - f24_im[v] * f35_im[v] );
                      par10_im[v] += - ( f24_re[v] * f35_im[v] + f24_im[v] * f35_re[v] );
                      // +[20][35]
                      par14_re[v] += + ( f20_re[v] * f35_re[v] - f20_im[v] * f35_im[v] );
                      par14_im[v] += + ( f20_re[v] * f35_im[v] + f20_im[v] * f35_re[v] );
                      // -[20][34]
                      par15_re[v] += - ( f20_re[v] * f34_re[v] - f20_im[v] * f34_im[v] );
                      par15_im[v] += - ( f20_re[v] * f34_im[v] + f20_im[v] * f34_re[v] );
                    }
                    
                    for(unsigned long int v=0; v<lv3; v++) {
                      corr.C[t][s0][s1][s2][s3][s4][s5][0] += sign * ( f10_re[v] * par10_re[v] - f10_im[v] * par10_im[v] );
                      corr.C[t][s0][s1][s2][s3][s4][s5][1] += sign * ( f10_re[v] * par10_im[v] + f10_im[v] * par10_re[v] );
                    }
                    for(unsigned long int v=0; v<lv3; v++) {
                      corr.C[t][s0][s1][s2][s3][s4][s5][0] += sign * ( f14_re[v] * par14_re[v] - f14_im[v] * par14_im[v] );
                      corr.C[t][s0][s1][s2][s3][s4][s5][1] += sign * ( f14_re[v] * par14_im[v] + f14_im[v] * par14_re[v] );
                    }
                    for(unsigned long int v=0; v<lv3; v++) {
                      corr.C[t][s0][s1][s2][s3][s4][s5][0] += sign * ( f15_re[v] * par15_re[v] - f15_im[v] * par15_im[v] );
                      corr.C[t][s0][s1][s2][s3][s4][s5][1] += sign * ( f15_re[v] * par15_im[v] + f15_im[v] * par15_re[v] );
                    }
                  }
                }
  }
}

inline void
qhg_fast_contract_f113(qhg_baryons_open_correlator corr, qhg_fast_spinor_field sp_1, qhg_fast_spinor_field sp_3) 
{

  qhg_lattice *lat = sp_1.lat; 
  int lt = lat->ldims[0];
  unsigned long int lv3 = lat->lv3;

#ifdef QHG_OMP
#pragma omp parallel for
#endif
  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++)
        for(int s2=0; s2<NS; s2++)
          for(int s3=0; s3<NS; s3++)
            for(int s4=0; s4<NS; s4++)
              for(int s5=0; s5<NS; s5++)
                for(int eps0=0; eps0<NEPS; eps0++) {
                  int c4 = QHG_EPS[eps0][0];
                  int c5 = QHG_EPS[eps0][1];
                  int c0 = QHG_EPS[eps0][2];
                  for(int eps1=0; eps1<NEPS; eps1++) {
                    int c2 = QHG_EPS[eps1][0];
                    int c3 = QHG_EPS[eps1][1];
                    int c1 = QHG_EPS[eps1][2];
                    double sign = QHG_EPS[eps0][3]*QHG_EPS[eps1][3];

                    /*
                     * NOTE: the loops have been packed in this way in order to induce self-vectorization
                     */

                    // ( [14][20] - [10][24] ) [35]
                    afloat par35_re[lv3];
                    afloat par35_im[lv3];

                    afloat * restrict f10_re = sp_1.field[t][s1][s0][c1][c0][0];
                    afloat * restrict f10_im = sp_1.field[t][s1][s0][c1][c0][1];
                    afloat * restrict f14_re = sp_1.field[t][s1][s4][c1][c4][0];
                    afloat * restrict f14_im = sp_1.field[t][s1][s4][c1][c4][1];
 
                    afloat * restrict f20_re = sp_1.field[t][s2][s0][c2][c0][0];
                    afloat * restrict f20_im = sp_1.field[t][s2][s0][c2][c0][1];
                    afloat * restrict f24_re = sp_1.field[t][s2][s4][c2][c4][0];
                    afloat * restrict f24_im = sp_1.field[t][s2][s4][c2][c4][1];
 
                    afloat * restrict f35_re = sp_3.field[t][s3][s5][c3][c5][0];
                    afloat * restrict f35_im = sp_3.field[t][s3][s5][c3][c5][1];

                    for(unsigned long int v=0; v<lv3; v++) {
                      // +[14][20]
                      par35_re[v]  = + ( f14_re[v] * f20_re[v] - f14_im[v] * f20_im[v] );
                      par35_im[v]  = + ( f14_re[v] * f20_im[v] + f14_im[v] * f20_re[v] );
                      
                      // -[10][24]
                      par35_re[v] += - ( f10_re[v] * f24_re[v] - f10_im[v] * f24_im[v] );
                      par35_im[v] += - ( f10_re[v] * f24_im[v] + f10_im[v] * f24_re[v] );
                    }

                    for(unsigned long int v=0; v<lv3; v++) {
                      corr.C[t][s0][s1][s2][s3][s4][s5][0] += sign * ( f35_re[v] * par35_re[v] - f35_im[v] * par35_im[v] );
                      corr.C[t][s0][s1][s2][s3][s4][s5][1] += sign * ( f35_re[v] * par35_im[v] + f35_im[v] * par35_re[v] );
                    }
                  }
                }
}

inline void
qhg_fast_contract_f121(qhg_baryons_open_correlator corr, qhg_fast_spinor_field sp_1, qhg_fast_spinor_field sp_2) 
{

  qhg_lattice *lat = sp_1.lat; 
  int lt = lat->ldims[0];
  unsigned long int lv3 = lat->lv3;

#ifdef QHG_OMP
#pragma omp parallel for
#endif
  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++)
        for(int s2=0; s2<NS; s2++)
          for(int s3=0; s3<NS; s3++)
            for(int s4=0; s4<NS; s4++)
              for(int s5=0; s5<NS; s5++)
                for(int eps0=0; eps0<NEPS; eps0++) {
                  int c4 = QHG_EPS[eps0][0];
                  int c5 = QHG_EPS[eps0][1];
                  int c0 = QHG_EPS[eps0][2];
                  for(int eps1=0; eps1<NEPS; eps1++) {
                    int c2 = QHG_EPS[eps1][0];
                    int c3 = QHG_EPS[eps1][1];
                    int c1 = QHG_EPS[eps1][2];
                    double sign = QHG_EPS[eps0][3]*QHG_EPS[eps1][3];

                    /*
                     * NOTE: the loops have been packed in this way in order to induce self-vectorization
                     */

                    // ( [15][30] - [10][35] ) [24]
                    afloat par24_re[lv3];
                    afloat par24_im[lv3];

                    afloat * restrict f10_re = sp_1.field[t][s1][s0][c1][c0][0];
                    afloat * restrict f10_im = sp_1.field[t][s1][s0][c1][c0][1];
                    afloat * restrict f15_re = sp_1.field[t][s1][s5][c1][c5][0];
                    afloat * restrict f15_im = sp_1.field[t][s1][s5][c1][c5][1];
 
                    afloat * restrict f30_re = sp_1.field[t][s3][s0][c3][c0][0];
                    afloat * restrict f30_im = sp_1.field[t][s3][s0][c3][c0][1];
                    afloat * restrict f35_re = sp_1.field[t][s3][s5][c3][c5][0];
                    afloat * restrict f35_im = sp_1.field[t][s3][s5][c3][c5][1];

                    afloat * restrict f24_re = sp_2.field[t][s2][s4][c2][c4][0];
                    afloat * restrict f24_im = sp_2.field[t][s2][s4][c2][c4][1];

                    for(unsigned long int v=0; v<lv3; v++) {
                      // +[15][30]
                      par24_re[v]  = + ( f15_re[v] * f30_re[v] - f15_im[v] * f30_im[v] );
                      par24_im[v]  = + ( f15_re[v] * f30_im[v] + f15_im[v] * f30_re[v] );
                      
                      // -[10][35]
                      par24_re[v] += - ( f10_re[v] * f35_re[v] - f10_im[v] * f35_im[v] );
                      par24_im[v] += - ( f10_re[v] * f35_im[v] + f10_im[v] * f35_re[v] );
                    }

                    for(unsigned long int v=0; v<lv3; v++) {
                      corr.C[t][s0][s1][s2][s3][s4][s5][0] += sign * ( f24_re[v] * par24_re[v] - f24_im[v] * par24_im[v] );
                      corr.C[t][s0][s1][s2][s3][s4][s5][1] += sign * ( f24_re[v] * par24_im[v] + f24_im[v] * par24_re[v] );
                    }
                  }
                }
}

inline void
qhg_fast_contract_f122(qhg_baryons_open_correlator corr, qhg_fast_spinor_field sp_1, qhg_fast_spinor_field sp_2) 
{

  qhg_lattice *lat = sp_1.lat; 
  int lt = lat->ldims[0];
  unsigned long int lv3 = lat->lv3;

#ifdef QHG_OMP
#pragma omp parallel for
#endif
  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++)
        for(int s2=0; s2<NS; s2++)
          for(int s3=0; s3<NS; s3++)
            for(int s4=0; s4<NS; s4++)
              for(int s5=0; s5<NS; s5++)
                for(int eps0=0; eps0<NEPS; eps0++) {
                  int c4 = QHG_EPS[eps0][0];
                  int c5 = QHG_EPS[eps0][1];
                  int c0 = QHG_EPS[eps0][2];
                  for(int eps1=0; eps1<NEPS; eps1++) {
                    int c2 = QHG_EPS[eps1][0];
                    int c3 = QHG_EPS[eps1][1];
                    int c1 = QHG_EPS[eps1][2];
                    double sign = QHG_EPS[eps0][3]*QHG_EPS[eps1][3];

                    /*
                     * NOTE: the loops have been packed in this way in order to induce self-vectorization
                     */
                    // [10] ( [25][34] - [24][35] )
                    afloat par10_re[lv3];
                    afloat par10_im[lv3];

                    // list of restrict pointers
                    afloat * restrict f10_re = sp_1.field[t][s1][s0][c1][c0][0];
                    afloat * restrict f10_im = sp_1.field[t][s1][s0][c1][c0][1];

                    afloat * restrict f24_re = sp_2.field[t][s2][s4][c2][c4][0];
                    afloat * restrict f24_im = sp_2.field[t][s2][s4][c2][c4][1];
                    afloat * restrict f25_re = sp_2.field[t][s2][s5][c2][c5][0];
                    afloat * restrict f25_im = sp_2.field[t][s2][s5][c2][c5][1];

                    afloat * restrict f34_re = sp_2.field[t][s3][s4][c3][c4][0];
                    afloat * restrict f34_im = sp_2.field[t][s3][s4][c3][c4][1];
                    afloat * restrict f35_re = sp_2.field[t][s3][s5][c3][c5][0];
                    afloat * restrict f35_im = sp_2.field[t][s3][s5][c3][c5][1];

                    for(unsigned long int v=0; v<lv3; v++) {
                      // +[25][34]
                      par10_re[v]  = + ( f25_re[v] * f34_re[v] - f25_im[v] * f34_im[v] );
                      par10_im[v]  = + ( f25_re[v] * f34_im[v] + f25_im[v] * f34_re[v] );

                      // -[24][35]
                      par10_re[v] += - ( f24_re[v] * f35_re[v] - f24_im[v] * f35_im[v] );
                      par10_im[v] += - ( f24_re[v] * f35_im[v] + f24_im[v] * f35_re[v] );

                    }

                    for(unsigned long int v=0; v<lv3; v++) {
                      corr.C[t][s0][s1][s2][s3][s4][s5][0] += sign * ( f10_re[v] * par10_re[v] - f10_im[v] * par10_im[v] );
                      corr.C[t][s0][s1][s2][s3][s4][s5][1] += sign * ( f10_re[v] * par10_im[v] + f10_im[v] * par10_re[v] );
                    }
                  }
                }
}

inline void
qhg_fast_contract_f123(qhg_baryons_open_correlator corr, qhg_fast_spinor_field sp_1, qhg_fast_spinor_field sp_2, qhg_fast_spinor_field sp_3) 
{

  qhg_lattice *lat = sp_1.lat; 
  int lt = lat->ldims[0];
  unsigned long int lv3 = lat->lv3;

#ifdef QHG_OMP
#pragma omp parallel for
#endif
  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++)
        for(int s2=0; s2<NS; s2++)
          for(int s3=0; s3<NS; s3++)
            for(int s4=0; s4<NS; s4++)
              for(int s5=0; s5<NS; s5++)
                for(int eps0=0; eps0<NEPS; eps0++) {
                  int c4 = QHG_EPS[eps0][0];
                  int c5 = QHG_EPS[eps0][1];
                  int c0 = QHG_EPS[eps0][2];
                  for(int eps1=0; eps1<NEPS; eps1++) {
                    int c2 = QHG_EPS[eps1][0];
                    int c3 = QHG_EPS[eps1][1];
                    int c1 = QHG_EPS[eps1][2];
                    double sign = QHG_EPS[eps0][3]*QHG_EPS[eps1][3];

                    /*
                     * NOTE: the loops have been packed in this way in order to induce self-vectorization
                     */

                    // - [10][24][35]
                    afloat par10_re[lv3];
                    afloat par10_im[lv3];

                    // list of restrict pointers
                    afloat * restrict f10_re = sp_1.field[t][s1][s0][c1][c0][0];
                    afloat * restrict f10_im = sp_1.field[t][s1][s0][c1][c0][1];

                    afloat * restrict f24_re = sp_2.field[t][s2][s4][c2][c4][0];
                    afloat * restrict f24_im = sp_2.field[t][s2][s4][c2][c4][1];

                    afloat * restrict f35_re = sp_3.field[t][s3][s5][c3][c5][0];
                    afloat * restrict f35_im = sp_3.field[t][s3][s5][c3][c5][1];

                    for(unsigned long int v=0; v<lv3; v++) {
                      // -[24][35]
                      par10_re[v]  = - ( f24_re[v] * f35_re[v] - f24_im[v] * f35_im[v] );
                      par10_im[v]  = - ( f24_re[v] * f35_im[v] + f24_im[v] * f35_re[v] );

                    }

                    for(unsigned long int v=0; v<lv3; v++) {
                      corr.C[t][s0][s1][s2][s3][s4][s5][0] += sign * ( f10_re[v] * par10_re[v] - f10_im[v] * par10_im[v] );
                      corr.C[t][s0][s1][s2][s3][s4][s5][1] += sign * ( f10_re[v] * par10_im[v] + f10_im[v] * par10_re[v] );
                    }
                  }
                }
}

qhg_baryons_open_correlator
qhg_fast_baryons(qhg_fast_spinor_field sp_1, qhg_fast_spinor_field sp_2, qhg_fast_spinor_field sp_3)
{
  qhg_lattice *lat = sp_1.lat; 
  int lt = lat->ldims[0];
  qhg_baryons_open_correlator corr;
  corr.C = qhg_alloc(lt*NS*NS*NS*NS*NS*NS*2*sizeof(double));
  memset(corr.C, '\0', lt*NS*NS*NS*NS*NS*NS*2*sizeof(double));
  corr.lat = lat; 

  // better to check "flavour" but not given at the moment
  // f111
  if(sp_1.field == sp_2.field && sp_1.field == sp_3.field)
    qhg_fast_contract_f111(corr, sp_1);
  // f113
  else if(sp_1.field == sp_2.field && sp_1.field != sp_3.field)
    qhg_fast_contract_f113(corr,sp_1,sp_3);
  // f121
  else if(sp_1.field != sp_2.field && sp_1.field == sp_3.field)
    qhg_fast_contract_f121(corr,sp_1,sp_2);
  // f122
  else if(sp_1.field != sp_2.field && sp_2.field == sp_3.field)
    qhg_fast_contract_f122(corr,sp_1,sp_2);
  // f123
  else if(sp_1.field != sp_2.field && sp_1.field != sp_3.field)
    qhg_fast_contract_f123(corr,sp_1,sp_2,sp_3);



  /* Split the communicator so that all processes of the same
     time-slices have a separate communicator. This way we can do a
     reduction over the new communicator. */
  MPI_Comm tcomm;
  int *pcoords = lat->comms->proc_coords;
  int proc_id = lat->comms->proc_id;
  MPI_Comm_split(lat->comms->comm, pcoords[0], proc_id, &tcomm);  

  /* The rank within the tcomm communicator */
  int s_proc_id;
  MPI_Comm_rank(tcomm, &s_proc_id);

  if(s_proc_id == 0) {
    MPI_Reduce( MPI_IN_PLACE, corr.C, lt*NS*NS*NS*NS*NS*NS*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    corr.write = 1;
  } else {
    MPI_Reduce( corr.C, NULL, lt*NS*NS*NS*NS*NS*NS*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    corr.write = 0;
  }

  //MPI_Comm_free(&tcomm);
  return corr;
}

void
qhg_baryons_open_correlator_finalize(qhg_baryons_open_correlator corr)
{
  free(corr.C);
  corr.lat = NULL;
  corr.write = 0;
}
  
void
qhg_write_baryons_open_correlator(char fname[], qhg_baryons_open_correlator corr, char group[])
{
  qhg_lattice *lat = corr.lat;
  int proc_id = lat->comms->proc_id;
  unsigned long int lvol = lat->lvol;
  int *pc = lat->comms->proc_coords;
  int *ld = lat->ldims;
  int *d = lat->dims;

  /* Split the communicator so that all processes that need to write are grouped together */
  MPI_Comm wcomm;
  MPI_Comm_split(lat->comms->comm, corr.write, proc_id, &wcomm);

  if ( corr.write == 0 )
    return;

  int ndims = 8;
  hsize_t starts[8] = {pc[0]*ld[0], 0, 0, 0, 0, 0, 0, 0};
  hsize_t dims[8] = {d[0], 4, 4, 4, 4, 4, 4, 2};
  hsize_t ldims[8] = {ld[0], 4, 4, 4, 4, 4, 4, 2};
  
  hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl_id, wcomm, MPI_INFO_NULL);
  hid_t file_id;
  if( access( fname, F_OK ) != -1 )
    file_id = H5Fopen(fname, H5F_ACC_RDWR, fapl_id);
  else
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
  H5Pclose(fapl_id);

  hid_t lcpl_id = H5Pcreate(H5P_LINK_CREATE);
  H5Pset_create_intermediate_group(lcpl_id, 1);  
  hid_t top_id = H5Gcreate(file_id, group, lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
  
  /*
    Attributes (metadata) are: 
    1) the origin (source position) and 
    2) the index order in the file
  */
  hsize_t n = 1;
  hid_t attrdat_id = H5Screate_simple(1, &n, NULL);
  hid_t attr_id = H5Acreate2(top_id, "Tsrc", H5T_NATIVE_INT, attrdat_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr_id, H5T_NATIVE_INT, &(corr.tsrc));
  H5Aclose(attr_id);
  H5Sclose(attrdat_id);

  char order[] = "C-order: [t,s0,s1,s3,s4,s5,s6,re/im]\0";
  attrdat_id = H5Screate(H5S_SCALAR);
  hid_t type_id = H5Tcopy(H5T_C_S1);
  H5Tset_size(type_id, strlen(order));
  attr_id = H5Acreate1(top_id, "IndexOrder", type_id, attrdat_id, H5P_DEFAULT);
  H5Awrite(attr_id, type_id, &order);

  H5Aclose(attr_id);
  H5Tclose(type_id);
  H5Sclose(attrdat_id);
  /* */
      
  hid_t filespace = H5Screate_simple(ndims, dims, NULL); 
  hid_t dataset_id = H5Dcreate(top_id, "corr", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t subspace = H5Screate_simple(ndims, ldims, NULL);
  filespace = H5Dget_space(dataset_id);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, starts, NULL, ldims, NULL);
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);      
  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, subspace, filespace, plist_id, corr.C);
  H5Dclose(dataset_id);
  H5Sclose(filespace);
  H5Sclose(subspace);
  H5Pclose(plist_id);
  H5Pclose(lcpl_id);
  H5Gclose(top_id);
  H5Fclose(file_id);
  MPI_Comm_free(&wcomm);
  return;
}
