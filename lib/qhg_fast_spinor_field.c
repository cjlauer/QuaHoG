#include <string.h>
#include <complex.h>
#include <qhg_alloc.h>
#include <qhg_defs.h>
#include <qhg_idx.h>
#include <qhg_types.h>
#include <qhg_correlator.h>
#include <qhg_xchange_utils.h>

qhg_fast_spinor_field
qhg_fast_spinor_field_init(qhg_lattice *lat, enum qhg_fermion_bc_time bc)
{
  unsigned long int lvol = lat->lvol;
  unsigned long int lv3 = lat->lv3;
  int lt = lat->ldims[0];
  
  qhg_fast_spinor_field sp;
  sp.alloc = qhg_alloc(lvol*NS*NS*NC*NC*2*sizeof(afloat));
  sp.field = qhg_alloc(lt*NS*NS*NC*NC*2*sizeof(afloat *));
  sp.lat = lat;
  sp.bc = bc;

  afloat *tmp = sp.alloc;

  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++)
        for(int c0=0; c0<NC; c0++)
          for(int c1=0; c1<NC; c1++)
            for(int ri=0; ri<2; ri++) {
              sp.field[t][s0][s1][c0][c1][ri] = tmp;
              tmp+=lv3;
            }

  return sp;
}

void
qhg_fast_spinor_field_finalize(qhg_fast_spinor_field sp)
{
  free(sp.alloc);
  free(sp.field);
  sp.lat = NULL;
  return;
}

void
qhg_fast_spinor_field_copy(qhg_fast_spinor_field y, qhg_fast_spinor_field x)
{
  y.lat = x.lat;
  y.bc = x.bc;
  memcpy(y.alloc, x.alloc, x.lat->lvol*NS*NS*NC*NC*2*sizeof(afloat));
  return;
}

void
qhg_fast_spinor_field_swap(qhg_fast_spinor_field *y, qhg_fast_spinor_field *x)
{
  y->lat = x->lat; 
  int lt = x->lat->ldims[0];
  enum qhg_fermion_bc_time bc = y->bc;
  y->bc = x->bc;
  x->bc = bc;
  afloat *tmp = y->alloc;
  y->alloc = x->alloc;
  x->alloc = tmp;
  
  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++)
        for(int c0=0; c0<NC; c0++)
          for(int c1=0; c1<NC; c1++)
            for(int ri=0; ri<2; ri++) {
              tmp = y->field[t][s0][s1][c0][c1][ri];
              y->field[t][s0][s1][c0][c1][ri] = x->field[t][s0][s1][c0][c1][ri];
              x->field[t][s0][s1][c0][c1][ri] = tmp;
            }
  return;
}

void
qhg_fast_spinor_field_import(qhg_fast_spinor_field y, qhg_spinor_field x[NS*NC])
{
  y.lat = x[0].lat;
  y.bc = x[0].bc;
  unsigned long int lvol = y.lat->lvol;
  unsigned long int lv3 = y.lat->lv3;

  for(int cs0=0; cs0<NS*NC; cs0++)
    for(unsigned long int v=0; v<lvol; v++)
      for(int cs1=0; cs1<NC*NS; cs1++) {
        y.field[v/lv3][CS2S(cs0)][CS2S(cs1)][CS2C(cs0)][CS2C(cs1)][0][v%lv3] = creal(x[cs1].field[(v*NC*NS + cs0)]);
        y.field[v/lv3][CS2S(cs0)][CS2S(cs1)][CS2C(cs0)][CS2C(cs1)][1][v%lv3] = cimag(x[cs1].field[(v*NC*NS + cs0)]);
      }
  return;
}

void
qhg_fast_spinor_field_import_propagator(qhg_fast_spinor_field y, qhg_spinor_field x, int cs1)
{
  y.lat = x.lat;
  y.bc = x.bc;
  unsigned long int lvol = y.lat->lvol;
  unsigned long int lv3 = y.lat->lv3;

  for(unsigned long int v=0; v<lvol; v++)
    for(int cs0=0; cs0<NC*NS; cs0++) {
      y.field[v/lv3][CS2S(cs0)][CS2S(cs1)][CS2C(cs0)][CS2C(cs1)][0][v%lv3] = creal(x.field[(v*NC*NS + cs0)]);
      y.field[v/lv3][CS2S(cs0)][CS2S(cs1)][CS2C(cs0)][CS2C(cs1)][1][v%lv3] = cimag(x.field[(v*NC*NS + cs0)]);
    }
  return;
}

void
qhg_fast_spinor_field_remove_bc_lower(qhg_fast_spinor_field x, int tsrc)
{
  switch(x.bc) {
  case PERIODIC:
    break;
  case ANTIPERIODIC:
    {
      int lt = x.lat->ldims[0];
      unsigned long int lv3 = x.lat->lv3;
      int t0 = x.lat->ldims[0] * x.lat->comms->proc_coords[0];
      unsigned long int lv = lv3*NC*NC*NS*NS*2;
      if ( t0 < tsrc ) {
        //printf("t0=%d src=%d max=%d\n", t0, tsrc, MIN(x.lat->ldims[0], tsrc-t0));
        for (int t = 0; t < (MIN(x.lat->ldims[0], tsrc-t0)); t++)
          for(int s0=0; s0<NS; s0++)
            for(int s1=0; s1<NS; s1++)
              for(int c0=0; c0<NC; c0++)
                for(int c1=0; c1<NC; c1++)
                  for(int ri=0; ri<2; ri++) {
                    afloat * restrict x_ = x.field[t][s0][s1][c0][c1][ri];
                    for(unsigned long int v=0; v<lv3; v++)
                      x_[v] *= -1.;
                  }
      }
      x.bc = PERIODIC;
    }
    break;
  }
  return;
}

void
qhg_fast_spinor_field_remove_bc_higher(qhg_fast_spinor_field x, int tsrc)
{
  switch(x.bc) {
  case PERIODIC:
    break;
  case ANTIPERIODIC:
    {
      int lt = x.lat->ldims[0];
      unsigned long int lv3 = x.lat->lv3;
      int t0 = x.lat->ldims[0] * x.lat->comms->proc_coords[0];
      unsigned long int lv = lv3*NC*NC*NS*NS*2;
      if ( t0 + lt > tsrc ) {
        //printf("t0=%d src=%d max=%d\n", t0, tsrc, MIN(x.lat->ldims[0], tsrc-t0));
        for (int t = (MAX(0,tsrc-t0)); t < lt; t++)
          for(int s0=0; s0<NS; s0++)
            for(int s1=0; s1<NS; s1++)
              for(int c0=0; c0<NC; c0++)
                for(int c1=0; c1<NC; c1++)
                  for(int ri=0; ri<2; ri++) {
                    afloat * restrict x_ = x.field[t][s0][s1][c0][c1][ri];
                    for(unsigned long int v=0; v<lv3; v++)
                      x_[v] *= -1.;
                  }
      }
      x.bc = PERIODIC;
    }
    break;
  }
  return;
}

#define SQUARE(a) (a)*(a)

double
qhg_fast_spinor_field_normsq(qhg_fast_spinor_field y)
{
  int lt = y.lat->ldims[0];
  unsigned long int lv3 = y.lat->lv3;

  double norm = 0;

  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++)
        for(int c0=0; c0<NC; c0++)
          for(int c1=0; c1<NC; c1++)
            for(int ri=0; ri<2; ri++) {
              afloat * restrict y_ = y.field[t][s0][s1][c0][c1][ri];
              for(unsigned long int v=0; v<lv3; v++) {
                norm += y_[v]*y_[v];
              }
            }

  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, y.lat->comms->comm);
  return norm;
}

void
qhg_fast_spinor_field_axpy(double a, qhg_fast_spinor_field x, qhg_fast_spinor_field y) {
  int lt = y.lat->ldims[0];
  unsigned long int lv3 = y.lat->lv3;
  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++)
        for(int c0=0; c0<NC; c0++)
          for(int c1=0; c1<NC; c1++) { 
            afloat * restrict x_re = x.field[t][s0][s1][c0][c1][0];
            afloat * restrict x_im = x.field[t][s0][s1][c0][c1][1];
            afloat * restrict y_re = y.field[t][s0][s1][c0][c1][0];
            afloat * restrict y_im = y.field[t][s0][s1][c0][c1][1];
            for(unsigned long int v=0; v<lv3; v++){
              x_re[v] = a*x_re[v] + y_re[v];
              x_im[v] = a*x_im[v] + y_im[v];
            }
          }
}

void
qhg_fast_spinor_field_trace_multiply_G_U_G(qhg_fast_spinor_field x, qhg_fast_gauge_field U, qhg_fast_spinor_field y, 
                                           qhg_correlator corr, int id, double factor)
{
  int lt = x.lat->ldims[0];
  unsigned long int lv3 = x.lat->lv3;
  unsigned long int lvol = x.lat->lvol;
  
  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++) {
        afloat z_re[lv3]; 
        afloat z_im[lv3]; 
        for(unsigned long int v=0; v<lv3; v++){
          z_re[v] = 0;
          z_im[v] = 0;
        }
        for(int s2=0; s2<NS; s2++)
          for(int c0=0; c0<NC; c0++)
            for(int c1=0; c1<NC; c1++) 
              for(int c2=0; c2<NC; c2++) {
                afloat * restrict x_re = x.field[t][s0][s2][c0][c1][0];
                afloat * restrict x_im = x.field[t][s0][s2][c0][c1][1];
                afloat * restrict U_re = U.field[t][c1][c2][0];
                afloat * restrict U_im = U.field[t][c1][c2][1];
                afloat * restrict y_re = y.field[t][s2][s1][c2][c0][0];
                afloat * restrict y_im = y.field[t][s2][s1][c2][c0][1];
                for(unsigned long int v=0; v<lv3; v++){
                  afloat re = x_re[v] * U_re[v] - x_im[v] * U_im[v];
                  afloat im = x_re[v] * U_im[v] + x_im[v] * U_re[v];
                  z_re[v] += y_re[v] * re - y_im[v] * im;
                  z_im[v] += y_re[v] * im + y_im[v] * re;
                }
            }
        for(unsigned long int v=0; v<lv3; v++){
          corr.C[(t*lv3 + v)*corr.site_size + (id*NS + s0)*NS +s1] = factor * (z_re[v] + I*z_im[v]);
        }
      }
}

void
qhg_fast_spinor_field_trace_multiply_U_G_G(qhg_fast_gauge_field U, qhg_fast_spinor_field x, qhg_fast_spinor_field y, 
                                           qhg_correlator corr, int id, double factor)
{
  int lt = x.lat->ldims[0];
  unsigned long int lv3 = x.lat->lv3;
  unsigned long int lvol = x.lat->lvol;
  
  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++) {
        afloat z_re[lv3]; 
        afloat z_im[lv3]; 
        for(unsigned long int v=0; v<lv3; v++){
          z_re[v] = 0;
          z_im[v] = 0;
        }
        for(int s2=0; s2<NS; s2++)
          for(int c0=0; c0<NC; c0++)
            for(int c1=0; c1<NC; c1++) 
              for(int c2=0; c2<NC; c2++) {
                afloat * restrict U_re = U.field[t][c0][c1][0];
                afloat * restrict U_im = U.field[t][c0][c1][1];
                afloat * restrict x_re = x.field[t][s0][s2][c1][c2][0];
                afloat * restrict x_im = x.field[t][s0][s2][c1][c2][1];
                afloat * restrict y_re = y.field[t][s2][s1][c2][c0][0];
                afloat * restrict y_im = y.field[t][s2][s1][c2][c0][1];
                for(unsigned long int v=0; v<lv3; v++){
                  afloat re = x_re[v] * U_re[v] - x_im[v] * U_im[v];
                  afloat im = x_re[v] * U_im[v] + x_im[v] * U_re[v];
                  z_re[v] += y_re[v] * re - y_im[v] * im;
                  z_im[v] += y_re[v] * im + y_im[v] * re;
                }
            }
        for(unsigned long int v=0; v<lv3; v++){
          corr.C[(t*lv3 + v)*corr.site_size + (id*NS + s0)*NS +s1] = factor * (z_re[v] + I*z_im[v]);
        }
      }
}

void
qhg_fast_spinor_field_multiply_G_G(qhg_fast_spinor_field z, qhg_fast_spinor_field x, qhg_fast_spinor_field y)
{
  int lt = x.lat->ldims[0];
  unsigned long int lv3 = x.lat->lv3;
  unsigned long int lvol = x.lat->lvol;
  
  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++) 
          for(int c0=0; c0<NC; c0++)
            for(int c1=0; c1<NC; c1++) { 
              afloat * restrict z_re = z.field[t][s0][s1][c0][c1][0];
              afloat * restrict z_im = z.field[t][s0][s1][c0][c1][1];
              for(unsigned long int v=0; v<lv3; v++){
                z_re[v] = 0;
                z_im[v] = 0;
              }
              for(int s2=0; s2<NS; s2++)
                for(int c2=0; c2<NC; c2++) {
                  afloat * restrict x_re = x.field[t][s0][s2][c0][c2][0];
                  afloat * restrict x_im = x.field[t][s0][s2][c0][c2][1];
                  afloat * restrict y_re = y.field[t][s2][s1][c2][c1][0];
                  afloat * restrict y_im = y.field[t][s2][s1][c2][c1][1];
                  for(unsigned long int v=0; v<lv3; v++){
                    z_re[v] += x_re[v] * y_re[v] - x_im[v] * y_im[v];
                    z_im[v] += x_re[v] * y_im[v] + x_im[v] * y_re[v];
                  }
                }
            }
}

void
qhg_fast_spinor_field_trace_multiply_U_G(qhg_fast_gauge_field U, qhg_fast_spinor_field x, 
                                           qhg_correlator corr, int id, double factor)
{
  int lt = x.lat->ldims[0];
  unsigned long int lv3 = x.lat->lv3;
  unsigned long int lvol = x.lat->lvol;
  
  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++) {
        afloat z_re[lv3]; 
        afloat z_im[lv3]; 
        for(unsigned long int v=0; v<lv3; v++){
          z_re[v] = 0;
          z_im[v] = 0;
        }
        for(int c0=0; c0<NC; c0++)
          for(int c1=0; c1<NC; c1++) {
            afloat * restrict U_re = U.field[t][c0][c1][0];
            afloat * restrict U_im = U.field[t][c0][c1][1];
            afloat * restrict x_re = x.field[t][s0][s1][c1][c0][0];
            afloat * restrict x_im = x.field[t][s0][s1][c1][c0][1];
            for(unsigned long int v=0; v<lv3; v++){
              z_re[v] += x_re[v] * U_re[v] - x_im[v] * U_im[v];
              z_im[v] += x_re[v] * U_im[v] + x_im[v] * U_re[v];
            }
          }
        for(unsigned long int v=0; v<lv3; v++){
          corr.C[(t*lv3 + v)*corr.site_size + (id*NS + s0)*NS +s1] = factor * (z_re[v] + I*z_im[v]);
        }
      }
}

void
qhg_fast_spinor_field_trace_multiply_Udag_G(qhg_fast_gauge_field U, qhg_fast_spinor_field x, 
                                           qhg_correlator corr, int id, double factor)
{
  int lt = x.lat->ldims[0];
  unsigned long int lv3 = x.lat->lv3;
  unsigned long int lvol = x.lat->lvol;
  
  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++) {
        afloat z_re[lv3]; 
        afloat z_im[lv3]; 
        for(unsigned long int v=0; v<lv3; v++){
          z_re[v] = 0;
          z_im[v] = 0;
        }
        for(int c0=0; c0<NC; c0++)
          for(int c1=0; c1<NC; c1++) { 
            afloat * restrict U_re = U.field[t][c1][c0][0];
            afloat * restrict U_im = U.field[t][c1][c0][1];
            afloat * restrict x_re = x.field[t][s0][s1][c1][c0][0];
            afloat * restrict x_im = x.field[t][s0][s1][c1][c0][1];
            for(unsigned long int v=0; v<lv3; v++){
              z_re[v] += x_re[v] * U_re[v] + x_im[v] * U_im[v];
              z_im[v] += x_im[v] * U_re[v] - x_re[v] * U_im[v];
            }
          }
        for(unsigned long int v=0; v<lv3; v++){
          corr.C[(t*lv3 + v)*corr.site_size + (id*NS + s0)*NS +s1] = factor * (z_re[v] + I*z_im[v]);
        }
      }
}

void
qhg_fast_spinor_field_trace_multiply_G_Udag_G(qhg_fast_spinor_field x, qhg_fast_gauge_field U, qhg_fast_spinor_field y,
                                              qhg_correlator corr, int id, double factor)
{
  int lt = x.lat->ldims[0];
  unsigned long int lv3 = x.lat->lv3;
  unsigned long int lvol = x.lat->lvol;

  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++) {
        afloat z_re[lv3]; 
        afloat z_im[lv3]; 
        for(unsigned long int v=0; v<lv3; v++){
          z_re[v] = 0;
          z_im[v] = 0;
        }
        for(int s2=0; s2<NS; s2++)
          for(int c0=0; c0<NC; c0++)
            for(int c1=0; c1<NC; c1++) 
              for(int c2=0; c2<NC; c2++) {
                afloat * restrict x_re = x.field[t][s0][s2][c0][c1][0];
                afloat * restrict x_im = x.field[t][s0][s2][c0][c1][1];
                afloat * restrict U_re = U.field[t][c2][c1][0];
                afloat * restrict U_im = U.field[t][c2][c1][1];
                afloat * restrict y_re = y.field[t][s2][s1][c2][c0][0];
                afloat * restrict y_im = y.field[t][s2][s1][c2][c0][1];
                for(unsigned long int v=0; v<lv3; v++){
                  afloat re = x_re[v] * U_re[v] + x_im[v] * U_im[v];
                  afloat im = - x_re[v] * U_im[v] + x_im[v] * U_re[v];
                  z_re[v] += y_re[v] * re - y_im[v] * im;
                  z_im[v] += y_re[v] * im + y_im[v] * re;
                }
            }
        for(unsigned long int v=0; v<lv3; v++){
          corr.C[(t*lv3 + v)*corr.site_size + (id*NS + s0)*NS +s1] = factor * (z_re[v] + I*z_im[v]);
        }
      }
}

void
qhg_fast_spinor_field_shift(qhg_fast_spinor_field y, qhg_fast_spinor_field x, int dir)
{
  y.lat = x.lat;
  y.bc = x.bc;
  int *ld = x.lat->ldims;
  int lt = x.lat->ldims[0];
  unsigned long int lv3 = x.lat->lv3;

  dir = (dir+ND)%(2*ND); // we change sign convention because the shift does the opposite thing compared to other shift functions.

  int fproc = x.lat->comms->neigh_proc[dir];
  int bproc = x.lat->comms->neigh_proc[(dir+ND)%(2*ND)];

  int sign = (dir / ND == 0) ? +1 : -1;
  dir = dir % ND;
  int proc_id = x.lat->comms->proc_id;
  MPI_Comm comm = x.lat->comms->comm;

  if (dir==0) {
    for(int t=0; t<lt; t++) {
      int t_mu = (t + lt + sign) % lt;
      for(int s0=0; s0<NS; s0++)
        for(int s1=0; s1<NS; s1++)
          for(int c0=0; c0<NC; c0++)
            for(int c1=0; c1<NC; c1++)
              for(int ri=0; ri<2; ri++) {
                afloat * restrict x_ = x.field[t][s0][s1][c0][c1][ri];
                afloat * restrict y_ = y.field[t_mu][s0][s1][c0][c1][ri];
                memcpy(y_, x_, lv3*sizeof(afloat));
              }
    }
    if(x.lat->comms->proc_dims[dir] > 1) {
      afloat *buffer;
      if(sign > 0)
        buffer = y.field[0][0][0][0][0][0];
      else
        buffer = y.field[lt-1][0][0][0][0][0];
      
      MPI_Sendrecv_replace(buffer, NS*NS*NC*NC*lv3*2*sizeof(afloat), MPI_BYTE, fproc, fproc, 
                           bproc, proc_id, comm, MPI_STATUS_IGNORE);
    }
  }
  else {
    unsigned long int v_out = lv3;
    for (int d = 1; d < dir; d++)
      v_out /= ld[d];
    unsigned long int v_in = v_out / ld[dir];
    afloat * buf = (x.lat->comms->proc_dims[dir] > 1) ? qhg_alloc(lt*lv3/ld[dir]*NS*NS*NC*NC*2*sizeof(afloat)) : NULL;
    afloat * restrict buf_ = buf;
    for(int t=0; t<lt; t++)
      for(int s0=0; s0<NS; s0++)
        for(int s1=0; s1<NS; s1++)
          for(int c0=0; c0<NC; c0++)
            for(int c1=0; c1<NC; c1++)
              for(int ri=0; ri<2; ri++) {
                afloat * restrict x_ = x.field[t][s0][s1][c0][c1][ri];
                afloat * restrict y_ = y.field[t][s0][s1][c0][c1][ri];
                if(x.lat->comms->proc_dims[dir] > 1) {
                  if(sign > 0) {
                    for(unsigned long int v=0; v<lv3/v_out; v++){
                      memcpy(y_+v_in, x_, (v_out-v_in)*sizeof(afloat));
                      memcpy(buf_, x_+(v_out-v_in), v_in*sizeof(afloat));
                      x_ += v_out;
                      y_ += v_out;
                      buf_ += v_in;
                    }
                  } else {
                    for(unsigned long int v=0; v<lv3/v_out; v++){
                      memcpy(y_, x_+v_in, (v_out-v_in)*sizeof(afloat));
                      memcpy(buf_, x_, v_in*sizeof(afloat));
                      x_ += v_out;
                      y_ += v_out;
                      buf_ += v_in;
                    }
                  }
                } else {
                  for(unsigned long int v=0; v<lv3/v_out; v++){
                    memcpy(y_+v_in, x_, (v_out-v_in)*sizeof(afloat));
                    memcpy(y_, x_+(v_out-v_in), v_in*sizeof(afloat));
                    x_ += v_out;
                    y_ += v_out;
                  }
                }
              }
    if(x.lat->comms->proc_dims[dir] > 1) {
      MPI_Sendrecv_replace(buf, lt*lv3/ld[dir]*NS*NS*NC*NC*2*sizeof(afloat), MPI_BYTE, fproc, fproc, 
                           bproc, proc_id, comm, MPI_STATUS_IGNORE);
      buf_ = buf;
      for(int t=0; t<lt; t++)
        for(int s0=0; s0<NS; s0++)
          for(int s1=0; s1<NS; s1++)
            for(int c0=0; c0<NC; c0++)
              for(int c1=0; c1<NC; c1++)
                for(int ri=0; ri<2; ri++) {
                  afloat * restrict y_ = y.field[t][s0][s1][c0][c1][ri];
                  if(sign > 0) {
                    for(unsigned long int v=0; v<lv3/v_out; v++){
                      memcpy(y_, buf_, v_in*sizeof(afloat));
                      y_ += v_out;
                      buf_ += v_in;
                    }
                  } else {
                    for(unsigned long int v=0; v<lv3/v_out; v++){
                      memcpy(y_+(v_out-v_in), buf_, v_in*sizeof(afloat));
                      y_ += v_out;
                      buf_ += v_in;
                    }
                  }
                }
      free(buf);
    }
  }
  return;
}
