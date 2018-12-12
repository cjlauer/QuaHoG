#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <complex.h>
#include <qhg_defs.h>
#include <qhg_idx.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_xchange_utils.h>
#include <qhg_correlator_shift.h>
/*
  This will swap every site with it's neighbor in +dir 

  1. If sign == +1, this will begin from the first site in direction
  "dir", and move forward in direction "dir". In the end the first
  slice will be the last
  
  2. If sign == -1, this will begin from the last slice in direction "dir"
  and move backwards in dir. In the end the last slice will be the first.  
 */
static void
shift_locally(_Complex double *c, size_t site_size, int dims[ND], int dir, int sign)
{
  _Complex double swap[site_size];  
  int d[] = {dims[0],dims[1],dims[2],dims[3]};
  d[dir] -= 1;
  
  for(int t=0; t<d[0]; t++)
    for(int x=0; x<d[1]; x++)
      for(int y=0; y<d[2]; y++)
	for(int z=0; z<d[3]; z++) {
	  long long int x0[] = {t,x,y,z};
	  long long int x1[] = {t,x,y,z};	    
	  if(sign == -1) {
	    x0[dir] = d[dir] - x0[dir] - 1;	    
	    x1[dir] = d[dir] - x1[dir] - 1;
	  }
	  x1[dir] += 1;
	  unsigned long int v1 = IDX(x1, dims);
	  unsigned long int v0 = IDX(x0, dims);
	  memcpy(swap, &c[v0*site_size], sizeof(_Complex double)*site_size);
	  memcpy(&c[v0*site_size], &c[v1*site_size], sizeof(_Complex double)*site_size);
	  memcpy(&c[v1*site_size], swap, sizeof(_Complex double)*site_size);	    
	}
  return;
}

void
qhg_correlator_shift(qhg_correlator corr, int shifts_in[ND])
{
  qhg_lattice *lat = corr.lat;
  int *dims = lat->dims;
  int shifts[ND];
  /* First ensure shift is within bounds */
  for(int dir=0; dir<ND; dir++)
    shifts[dir] = (shifts_in[dir] + dims[dir]) % dims[dir];
  
  /* Get the shortest absolute distance */
  for(int dir=0; dir<ND; dir++)
    shifts[dir] = shifts[dir] <= dims[dir]/2 ? shifts[dir] : (shifts[dir] - dims[dir]);
  
  size_t site_size = corr.site_size;
  int *ldims = lat->ldims;
  int *procs = lat->comms->proc_dims;
  MPI_Comm comm = lat->comms->comm;
  int proc_id = lat->comms->proc_id;  
  size_t site_size_bytes = site_size*sizeof(_Complex double);
  for(int dir=0; dir<ND; dir++) {    
    int sign = shifts[dir] >= 0 ? +1 : -1;
    MPI_Datatype send_type = get_slice(ND, ldims, dir, sign == +1 ? 0 : ldims[dir]-1, site_size_bytes);
    int fproc = lat->comms->neigh_proc[dir];
    int bproc = lat->comms->neigh_proc[dir+ND];
    for(int i=0; i<abs(shifts[dir]); i++) {
      if(procs[dir] > 1) {
	if(sign > 0) {
	  /* 
	     Replace boundary face dir = 0 with boundary face dir = 0
	     from forward neighbor
	   */
	  MPI_Sendrecv_replace(corr.C, 1, send_type, bproc, bproc, fproc, proc_id,
			       comm, MPI_STATUS_IGNORE);
	} else {
	  /* 
	     Replace boundary face dir = L-1 with boundary face dir =
	     L-1 from backwards neighbor
	   */
	  MPI_Sendrecv_replace(corr.C, 1, send_type, fproc, fproc, bproc, proc_id,
			       comm, MPI_STATUS_IGNORE);
	}
      }      
      shift_locally(corr.C, site_size, ldims, dir, sign);
    }
    MPI_Type_free(&send_type);
  }

  /*
    Correct origin by shift
  */
  for(int dir=0; dir<ND; dir++)
    corr.origin[dir] = (dims[dir] + (corr.origin[dir] - shifts[dir])) % dims[dir];

  return;
}

static void
shift_locally_ft(_Complex double *c, size_t site_size, int lt, int sign)
{
  _Complex double swap[site_size];  
  
  for(int t=0; t<lt-1; t++) {
    int t0 = t;
    int t1 = t;	    
    if(sign == -1) {
      t0 = lt - t0 - 2;	    
      t1 = lt - t1 - 2;	    
    }
    t1 += 1;
    memcpy(swap, &c[t0*site_size], sizeof(_Complex double)*site_size);
    memcpy(&c[t0*site_size], &c[t1*site_size], sizeof(_Complex double)*site_size);
    memcpy(&c[t1*site_size], swap, sizeof(_Complex double)*site_size);	    
  }
  return;
}

void
qhg_correlator_shift_ft(qhg_correlator corr, int shift, MPI_Comm time_comm)
{
  if(shift==0) return;
  qhg_lattice *lat = corr.lat;
  int lt = lat->dims[0];
  /* First ensure shift is within bounds */
  shift = (shift + lt) % lt;
  /* Get the shortest absolute distance */
  shift = shift <= lt/2 ? shift : (shift - lt);
  
  size_t site_size = corr.site_size * corr.mom_list->n_mom_vecs;

  int proc_id = lat->comms->proc_id;  
  size_t site_size_bytes = site_size*sizeof(_Complex double);
  int sign = shift >= 0 ? +1 : -1;
  MPI_Datatype send_type = get_slice(1, &lt, 0, sign == +1 ? 0 : lt-1, site_size_bytes);
  int time_rank, time_size;
  MPI_Comm_size( time_comm, &time_size);
  MPI_Comm_rank( time_comm, &time_rank );
  int fproc = (time_rank + 1)%time_size;
  int bproc = (time_rank + time_size - 1)%time_size;
  for(int i=0; i<abs(shift); i++) {
    if(time_size > 1) {
      if(sign > 0) {
        /* 
           Replace boundary face dir = 0 with boundary face dir = 0
           from forward neighbor
        */
        MPI_Sendrecv_replace(corr.C, 1, send_type, bproc, bproc, fproc, time_rank,
                             time_comm, MPI_STATUS_IGNORE);
      } else {
        /* 
           Replace boundary face dir = L-1 with boundary face dir =
           L-1 from backwards neighbor
        */
        MPI_Sendrecv_replace(corr.C, 1, send_type, fproc, fproc, bproc, time_rank,
                             time_comm, MPI_STATUS_IGNORE);
      }
    }
    shift_locally_ft(corr.C, site_size, lt, sign);
  }
  MPI_Type_free(&send_type);

  /*
    Correct origin by shift
  */
  corr.origin[0] = (lt + (corr.origin[0] - shift)) % lt;

  return;
}
