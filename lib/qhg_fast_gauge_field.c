#include <mpi.h>
#include <string.h>
#include <complex.h>
#include <qhg_alloc.h>
#include <qhg_defs.h>
#include <qhg_idx.h>
#include <qhg_types.h>
#include <qhg_xchange_utils.h>

qhg_fast_gauge_field
qhg_fast_gauge_field_init(qhg_lattice *lat)
{
  unsigned long int lvol = lat->lvol;
  unsigned long int lv3 = lat->lv3;
  int lt = lat->ldims[0];
  
  qhg_fast_gauge_field gf;
  gf.alloc = qhg_alloc(lvol*NC*NC*2*sizeof(afloat));
  gf.field = qhg_alloc(lt*NC*NC*2*sizeof(afloat *));
  gf.lat = lat;
  
  afloat *tmp = gf.alloc;
  
  for(int t=0; t<lt; t++)
    for(int c0=0; c0<NC; c0++)
      for(int c1=0; c1<NC; c1++)
        for(int ri=0; ri<2; ri++) {
          gf.field[t][c0][c1][ri] = tmp;
          tmp+=lv3;
        }

  return gf;
}

void
qhg_fast_gauge_field_finalize(qhg_fast_gauge_field gf)
{
  free(gf.alloc);
  free(gf.field);
  gf.lat = NULL;
  return;
}

void
qhg_fast_gauge_field_copy(qhg_fast_gauge_field y, qhg_fast_gauge_field x)
{
  y.lat = x.lat;
  memcpy(y.alloc, x.alloc, x.lat->lvol*NC*NC*2*sizeof(afloat));
  return;
}

// The link along dir is imported exiting the site
void
qhg_fast_gauge_field_import_dir(qhg_fast_gauge_field y, qhg_gauge_field x, int dir)
{
  int sign = dir / ND == 0 ? +1 : -1;
  dir = dir % ND;

  y.lat = x.lat;
  int lt = x.lat->ldims[0];
  unsigned long int lv3 = x.lat->lv3;

  if ( sign < 0 ) {
    // we copy from the neighbor the link in direction dir, daggered
    _Complex double * xf = x.field + dir*NC*NC;
    for(int t=0; t<lt; t++)
      for(unsigned long int v=0; v<lv3; v++){
        unsigned long int v_mu = x.lat->nn[dir+ND][t*lv3+v];
        _Complex double * xf = x.field + v_mu*NC*NC*ND + dir*NC*NC;
        for(int c0=0; c0<NC; c0++)
          for(int c1=0; c1<NC; c1++) {
            y.field[t][c1][c0][0][v] = creal(*xf);
            y.field[t][c1][c0][1][v] = -cimag(*xf);
            xf++;
          }
      }
  }
  else {
    // we copy in-place the link in direction dir
    _Complex double * xf = x.field + dir*NC*NC;
    for(int t=0; t<lt; t++)
      for(unsigned long int v=0; v<lv3; v++){
        for(int c0=0; c0<NC; c0++)
          for(int c1=0; c1<NC; c1++) {
            y.field[t][c0][c1][0][v] = creal(*xf);
            y.field[t][c0][c1][1][v] = cimag(*xf);
            xf++;
          }
        xf += (ND-1)*NC*NC;
      }
  }
  return;
}

// The link along dir is imported entering the site
void
qhg_fast_gauge_field_import_dir_daggered(qhg_fast_gauge_field y, qhg_gauge_field x, int dir)
{
  int sign = dir / ND == 0 ? +1 : -1;
  dir = dir % ND;

  y.lat = x.lat;
  int lt = x.lat->ldims[0];
  unsigned long int lv3 = x.lat->lv3;

  if ( sign < 0 ) {
    // we copy from the neighbor the link in direction dir
    for(int t=0; t<lt; t++)
      for(unsigned long int v=0; v<lv3; v++){
        unsigned long int v_mu = x.lat->nn[dir+ND][t*lv3+v];
        _Complex double * xf = x.field + v_mu*NC*NC*ND + dir*NC*NC;
        for(int c0=0; c0<NC; c0++)
          for(int c1=0; c1<NC; c1++) {
            y.field[t][c0][c1][0][v] = creal(*xf);
            y.field[t][c0][c1][1][v] = cimag(*xf);
            xf++;
          }
      }
  }
  else {
    // we copy in-place the link in direction dir, daggered
    _Complex double * xf = x.field + dir*NC*NC;
    for(int t=0; t<lt; t++)
      for(unsigned long int v=0; v<lv3; v++){
        for(int c0=0; c0<NC; c0++)
          for(int c1=0; c1<NC; c1++) {
            y.field[t][c1][c0][0][v] = creal(*xf);
            y.field[t][c1][c0][1][v] = -cimag(*xf);
            xf++;
          }
        xf += (ND-1)*NC*NC;
      }
  }
  return;
}

void
qhg_fast_gauge_field_multiply(qhg_fast_gauge_field z, qhg_fast_gauge_field x, qhg_fast_gauge_field y)
{
  z.lat = x.lat;
  int lt = x.lat->ldims[0];
  unsigned long int lv3 = x.lat->lv3;
  unsigned long int lvol = x.lat->lvol;

  memset(z.alloc, '\0', lvol*NC*NC*2*sizeof(afloat));

  for(int t=0; t<lt; t++)
    for(int c0=0; c0<NC; c0++)
      for(int c1=0; c1<NC; c1++) 
        for(int c2=0; c2<NC; c2++) {
          afloat * restrict x_re = x.field[t][c0][c1][0];
          afloat * restrict x_im = x.field[t][c0][c1][1];
          afloat * restrict y_re = y.field[t][c1][c2][0];
          afloat * restrict y_im = y.field[t][c1][c2][1];
          afloat * restrict z_re = z.field[t][c0][c2][0];
          afloat * restrict z_im = z.field[t][c0][c2][1];
          for(unsigned long int v=0; v<lv3; v++){
            z_re[v] += x_re[v] * y_re[v] - x_im[v] * y_im[v];
            z_im[v] += x_re[v] * y_im[v] + x_im[v] * y_re[v];
          }
        }
  return;
}

void
qhg_fast_gauge_field_id(qhg_fast_gauge_field x) {
  unsigned long int lv3 = x.lat->lv3;
  int lt = x.lat->ldims[0];
  for(int t=0; t<lt; t++)
    for(int c0=0; c0<NC; c0++)
      for(int c1=0; c1<NC; c1++) {
        afloat * restrict x_re = x.field[t][c0][c1][0];
        afloat * restrict x_im = x.field[t][c0][c1][1];
        if (c0==c1)
          for(unsigned long int v=0; v<lv3; v++)
            x_re[v] = 1;
        else
          memset(x_re, '\0', lv3*sizeof(afloat));
        memset(x_im, '\0', lv3*sizeof(afloat));
      }
}

void
qhg_fast_gauge_field_shift(qhg_fast_gauge_field y, qhg_fast_gauge_field x, int dir)
{
  y.lat = x.lat;
  int *ld = x.lat->ldims;
  int lt = x.lat->ldims[0];
  unsigned long int lv3 = x.lat->lv3;

  dir = (dir+ND)%(2*ND); // we change sign convention because the shift does the opposite thing compared to other shift functions.

  int fproc = x.lat->comms->neigh_proc[dir];
  int bproc = x.lat->comms->neigh_proc[(dir+ND)%(2*ND)];

  int sign = dir / ND == 0 ? +1 : -1;
  dir = dir % ND;
  int proc_id = x.lat->comms->proc_id;
  MPI_Comm comm = x.lat->comms->comm;

  if (dir==0) {
    for(int t=0; t<lt; t++) {
      int t_mu = (t + lt + sign) % lt;
      for(int c0=0; c0<NC; c0++)
        for(int c1=0; c1<NC; c1++)
          for(int ri=0; ri<2; ri++) {
            afloat * restrict x_ = x.field[t][c0][c1][ri];
            afloat * restrict y_ = y.field[t_mu][c0][c1][ri];
            memcpy(y_, x_, lv3*sizeof(afloat));
          }
    }
    if(x.lat->comms->proc_dims[dir] > 1) {
      afloat *buffer;
      if(sign > 0)
        buffer = y.field[0][0][0][0];
      else
        buffer = y.field[lt-1][0][0][0];
      
      MPI_Sendrecv_replace(buffer, NC*NC*lv3*2*sizeof(afloat), MPI_BYTE, fproc, fproc, 
                           bproc, proc_id, comm, MPI_STATUS_IGNORE);
    }
  }
  else {
    unsigned long int v_out = lv3;
    for (int d = 1; d < dir; d++)
      v_out /= ld[d];
    unsigned long int v_in = v_out / ld[dir];
    afloat * buf = (x.lat->comms->proc_dims[dir] > 1) ? qhg_alloc(lt*lv3/ld[dir]*NC*NC*2*sizeof(afloat)) : NULL;
    afloat * restrict buf_ = buf;
    for(int t=0; t<lt; t++)
      for(int c0=0; c0<NC; c0++)
        for(int c1=0; c1<NC; c1++)
          for(int ri=0; ri<2; ri++) {
            afloat * restrict x_ = x.field[t][c0][c1][ri];
            afloat * restrict y_ = y.field[t][c0][c1][ri];
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
      MPI_Sendrecv_replace(buf, lt*lv3/ld[dir]*NC*NC*2*sizeof(afloat), MPI_BYTE, fproc, fproc, 
                           bproc, proc_id, comm, MPI_STATUS_IGNORE);
      buf_ = buf;
      
      for(int t=0; t<lt; t++)
        for(int c0=0; c0<NC; c0++)
          for(int c1=0; c1<NC; c1++)
            for(int ri=0; ri<2; ri++) {
              afloat * restrict y_ = y.field[t][c0][c1][ri];
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
