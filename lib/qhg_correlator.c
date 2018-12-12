#include <string.h>

#include <qhg_alloc.h>
#include <qhg_defs.h>
#include <qhg_types.h>

qhg_correlator
qhg_correlator_init(size_t site_size, qhg_lattice *lat)
{
  unsigned long int lvol = lat->lvol;  
  qhg_correlator corr;
  corr.C = qhg_alloc(lvol*site_size*sizeof(_Complex double));
  memset(corr.C, '\0', lvol*site_size*sizeof(_Complex double));
  corr.lat = lat;
  corr.mom_list = NULL;
  corr.site_size = site_size;
  corr.origin = qhg_alloc(ND*sizeof(int));
  for(int i=0; i<ND; i++)
    corr.origin[i] = 0;
  return corr;
}

qhg_correlator
qhg_transformed_correlator_init(size_t site_size, qhg_lattice *lat, qhg_mom_list *mom_list)
{
  unsigned long int n_mom = mom_list->n_mom_vecs;  
  qhg_correlator corr;
  corr.C = qhg_alloc(lat->dims[0]*n_mom*site_size*sizeof(_Complex double));
  memset(corr.C, '\0', lat->dims[0]*n_mom*site_size*sizeof(_Complex double));
  corr.lat = lat;
  corr.mom_list = mom_list;
  corr.site_size = site_size;
  corr.origin = qhg_alloc(1*sizeof(int));
  for(int i=0; i<1; i++)
    corr.origin[i] = 0;
  return corr;
}

qhg_correlator
qhg_correlator_copy(qhg_correlator x)
{
  qhg_lattice *lat = x.lat;
  size_t site_size = x.site_size;
  qhg_correlator y;
  unsigned long int n_mom = x.mom_list->n_mom_vecs;

  if (x.mom_list == NULL) {
    y = qhg_correlator_init(site_size, lat);  
    memcpy(y.C, x.C, lat->lvol*site_size*sizeof(_Complex double));
    for(int i=0; i<ND; i++)
      y.origin[i] = x.origin[i];
  }
  else {
    y = qhg_transformed_correlator_init(site_size, lat, x.mom_list);  
    memcpy(y.C, x.C, lat->dims[0]*n_mom*site_size*sizeof(_Complex double));
    for(int i=0; i<1; i++)
      y.origin[i] = x.origin[i];
  }
  return y;
}

void
qhg_correlator_add(qhg_correlator x, qhg_correlator y)
{
  size_t vol = x.lat->lvol;
  size_t site_size = x.site_size;

  for(size_t i=0; i<vol*site_size; i++) {
    x.C[i] += y.C[i];
  }
}

void
qhg_correlator_finalize(qhg_correlator corr)
{
  free(corr.C);
  corr.lat = NULL;
  free(corr.origin);
  corr.origin = NULL;
  return;
}
