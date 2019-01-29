#include <string.h>
#include <stdio.h>
#include <qhg_alloc.h>
#include <qhg_defs.h>
#include <qhg_types.h>

qhg_der_correlator
qhg_der_correlator_init(size_t site_size, qhg_lattice *lat, qhg_mom_list *mom_list, int der_order)
{

  qhg_der_correlator corr;
  int nmoms;
  int (*moms)[3];
  corr.lat = lat;
  corr.der_order = der_order;
  corr.site_size = site_size;
  corr.ncorr = 1<<(2*der_order);
  if(mom_list == NULL) {
    corr.vol_size = lat->lvol;
    corr.mom_list = NULL;
    corr.origin = qhg_alloc(ND*sizeof(int));
    for(int i=0; i<ND; i++)
      corr.origin[i] = 0;
  } else {
    corr.mom_list = mom_list;
    nmoms = corr.mom_list->n_mom_vecs;
    moms = corr.mom_list->mom_vecs;
    corr.vol_size = lat->dims[0]*nmoms;
    corr.origin = qhg_alloc(sizeof(int));
    corr.origin = 0;
  }
  corr.C = qhg_alloc(corr.ncorr*sizeof(_Complex double *));
  for(int i=0; i<corr.ncorr; i++) { 
    corr.C[i] = qhg_alloc(corr.vol_size*site_size*sizeof(_Complex double));
    memset(corr.C[i], '\0', corr.vol_size*site_size*sizeof(_Complex double));
  }

  return corr;
}

qhg_der_correlator
qhg_averaged_der_correlator_init(size_t site_size, qhg_lattice *lat, qhg_mom_list *mom_list, int ncorr, int der_order)
{

  qhg_der_correlator corr;
  int nmoms;
  int (*moms)[3];
  corr.lat = lat;
  corr.der_order = der_order;
  corr.site_size = site_size;
  corr.ncorr = ncorr;
  if(mom_list == NULL) {
    corr.vol_size = lat->lvol;
    corr.mom_list = NULL;
    corr.origin = qhg_alloc(ND*sizeof(int));
    for(int i=0; i<ND; i++)
      corr.origin[i] = 0;
  } else {
    corr.mom_list = mom_list;
    nmoms = corr.mom_list->n_mom_vecs;
    moms = corr.mom_list->mom_vecs;
    corr.vol_size = lat->dims[0]*nmoms;
    corr.origin = qhg_alloc(sizeof(int));
    corr.origin = 0;
  }
  corr.C = qhg_alloc(ncorr*sizeof(_Complex double *));
  for(int i=0; i<ncorr; i++) { 
    corr.C[i] = qhg_alloc(corr.vol_size*site_size*sizeof(_Complex double));
    memset(corr.C[i], '\0', corr.vol_size*site_size*sizeof(_Complex double));
  }

  return corr;
}

qhg_der_correlator
qhg_der_correlator_copy(qhg_der_correlator x)
{
  qhg_lattice *lat = x.lat;
  size_t site_size = x.site_size;
  qhg_der_correlator y;
  qhg_mom_list *mom_list = x.mom_list;
  int der_order = x.der_order;
  int ncorr = x.ncorr;

  y = qhg_der_correlator_init(site_size, lat, mom_list, der_order);  

  y.dt = x.dt;
  y.der_order = der_order;
  y.ncorr = ncorr;
  y.flav = x.flav;
  y.proj = x.proj;

  for(int i=0; i<ncorr; i++) {
    memcpy(y.C[i], x.C[i], x.vol_size*site_size*sizeof(_Complex double));
  }
  if (mom_list == NULL) {
    for(int i=0; i<ND; i++)
      y.origin[i] = x.origin[i];
  }
  else {
    y.origin[0] = x.origin[0];
  }

  return y;
}

qhg_der_correlator
qhg_averaged_der_correlator_copy(qhg_der_correlator x)
{
  qhg_lattice *lat = x.lat;
  size_t site_size = x.site_size;
  qhg_der_correlator y;
  qhg_mom_list *mom_list = x.mom_list;
  int der_order = x.der_order;
  int ncorr = x.ncorr;

  y = qhg_averaged_der_correlator_init(site_size, lat, mom_list, ncorr, der_order);  

  y.dt = x.dt;
  y.der_order = der_order;
  y.ncorr = ncorr;
  y.flav = x.flav;
  y.proj = x.proj;

  for(int i=0; i<ncorr; i++) {
    memcpy(y.C[i], x.C[i], x.vol_size*site_size*sizeof(_Complex double));
  }
  if (mom_list == NULL) {
    for(int i=0; i<ND; i++)
      y.origin[i] = x.origin[i];
  }
  else {
    y.origin[0] = x.origin[0];
  }

  return y;
}

void
qhg_der_correlator_finalize(qhg_der_correlator corr)
{
  for(int i=0; i<corr.ncorr; i++)
    free(corr.C[i]);
    corr.C = NULL;
  corr.lat = NULL;
  free(corr.origin);
  corr.origin = NULL;
  return;
}
