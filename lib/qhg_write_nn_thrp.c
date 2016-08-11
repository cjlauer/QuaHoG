#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

#include <complex.h>
#include <qhg_idx.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_correlator.h>
#include <qhg_prop_gammas.h>
#include <qhg_prop_ops.h>
#include <qhg_io_utils.h>
#include <qhg_nn_thrp_defs.h>

void
qhg_write_nn_thrp(char fname[], qhg_thrp_correlator corr_thrp, char group[])
{
  qhg_lattice *lat = corr_thrp.corr.lat;
  int proc_id = lat->comms->proc_id;
  int am_io_proc = proc_id == 0 ? 1 : 0;
  int site_size = corr_thrp.corr.site_size;
  int *pd = lat->comms->proc_dims;
  int *pc = lat->comms->proc_coords;
  int np = lat->comms->nprocs;
  int *ld = lat->ldims;
  int *d = lat->dims;

  int offset = 0;
  hsize_t starts[ND+1] = {pc[0]*ld[0], pc[1]*ld[1], pc[2]*ld[2], pc[3]*ld[3], 0};
  hsize_t dims[ND+1] = {d[0], d[1], d[2], d[3], 2};
  hsize_t ldims[ND+1] = {ld[0], ld[1], ld[2], ld[3], 2};

  /*
    We will only save the set of time-slices between [source, sink]
  */
  {
    int t0 = corr_thrp.corr.origin[0];
    int dt = corr_thrp.dt;
    int ts = t0 + dt;
    int Lt = dims[0];
    int lt = ldims[0];
    int t_proc_t0 = t0/lt;		/* t-rank of processes with the source time-slice */
    int t_proc_ts = ((ts%Lt)/lt)%pd[0]; /* t-rank of processes with the sink time-slice */
    /* dims[] and ldims[] will be the global and local dimensions of the source-to-sink data chunck */
    dims[0] = dt+1;
    ldims[0] = 0;
    if(pc[0] == t_proc_t0)
      ldims[0] = lt;

    /* If this is one of the processes with the sink time-slice */
    if(pc[0] == t_proc_ts)
      ldims[0] = lt;
    
    /* Start from the ranks with t0 and iterate up to the ranks with ts */
    int r = t_proc_t0;
    while(r != t_proc_ts) {
      if(pc[0] == r)
  	ldims[0] = lt;
      r = (r+1) % pd[0];
    }
    /* All participating ranks now have the local t-dimension set to
       lt, for everyone else it is set to 0. Now adjust the local
       t-dimension of the edge ranks, i.e. the one with t0 and the one
       with ts */
    if(pc[0] == t_proc_t0)
      ldims[0] = lt-t0%lt;
    
    if(pc[0] == t_proc_ts)
      ldims[0] = ldims[0] - (lt-ts%lt) + 1;

    if((pc[0] == t_proc_ts) && (pc[0] == t_proc_t0))
      ldims[0] = dt+1;
    
    /* Now find the t-index of the chunck that is to be written which
       each rank starts from */
    starts[0] = (Lt + lt*pc[0]-t0) % Lt;

    /* This is now correct for all participating ranks except for
       those which hold t0, for which the start should be zero, and we
       need to take care in shifting the data to the origin when
       copying to buf[] further down */
    if(pc[0] == t_proc_t0) {
      offset = t0 % lt;
      starts[0] = 0;
    }
  }
  
  hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
  H5Pclose(fapl_id);

  size_t lvol = 1;
  for(int i=0; i<ND; i++)
    lvol *= ldims[i];

  hid_t lcpl_id = H5Pcreate(H5P_LINK_CREATE);
  H5Pset_create_intermediate_group(lcpl_id, 1);  
  hid_t top_id = H5Gcreate(file_id, group, lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
  
  double *buf = qhg_alloc(sizeof(double)*lvol*2);
  _Complex double *c = corr_thrp.corr.C;
  /*
    Hierarchy is /insertion/ with insertions listed in the
    qhg_nn_thrp_defs.h header file
  */
  for(int ichan=0; ichan<NCHAN; ichan++) {
    for(int t0=0; t0<ldims[0]; t0++) {
      int lt = lat->ldims[0];
      int lv3 = lat->lvol/lt;
      int t = (t0+offset)%lt;
      for(int v3=0; v3<lv3; v3++) {
	int v0 = t0*lv3 + v3;
	int v = t*lv3 + v3;
	buf[0 + 2*v0] = creal(c[TIDX(v, ichan)]);
	buf[1 + 2*v0] = cimag(c[TIDX(v, ichan)]);
      }
    }

    char *group_tag;
    asprintf(&group_tag, "%s", chan_tags[ichan]);
    hid_t group_id = H5Gcreate(top_id, group_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
      Attributes (metadata) are:
      1) the origin (source position) [1-dimensional integer array of 4 elements]
      2) the index order in the file [string]
      3) source-sink separation [single integer]
      4) flavor [string]
      5) projector [string]
    */

    {
      hsize_t n = 4;
      hid_t attrdat_id = H5Screate_simple(1, &n, NULL);
      hid_t attr_id = H5Acreate2(group_id, "Origin", H5T_NATIVE_INT, attrdat_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attr_id, H5T_NATIVE_INT, corr_thrp.corr.origin);
      H5Aclose(attr_id);
      H5Sclose(attrdat_id);
    }

    {
      char order[] = "C-order: [t,x,y,z,real/imag]\0";
      hid_t attrdat_id = H5Screate(H5S_SCALAR);
      hid_t type_id = H5Tcopy(H5T_C_S1);
      H5Tset_size(type_id, strlen(order));
      hid_t attr_id = H5Acreate1(group_id, "IndexOrder", type_id, attrdat_id, H5P_DEFAULT);
      H5Awrite(attr_id, type_id, &order);
      H5Aclose(attr_id);
      H5Tclose(type_id);
      H5Sclose(attrdat_id);
    }

    {
      hid_t attrdat_id = H5Screate(H5S_SCALAR);
      hid_t attr_id = H5Acreate1(group_id, "Sink-source separation", H5T_NATIVE_INT, attrdat_id, H5P_DEFAULT);
      H5Awrite(attr_id, H5T_NATIVE_INT, &corr_thrp.dt);
      H5Aclose(attr_id);
      H5Sclose(attrdat_id);
    }
  
    {
      char *flav;
      switch(corr_thrp.flav) {
      case up:
	flav = strdup("up");
	break;
      case dn:
	flav = strdup("dn");
	break;
      }      
      hid_t attrdat_id = H5Screate(H5S_SCALAR);
      hid_t type_id = H5Tcopy(H5T_C_S1);
      H5Tset_size(type_id, strlen(flav));
      hid_t attr_id = H5Acreate1(group_id, "Flavor", type_id, attrdat_id, H5P_DEFAULT);
      H5Awrite(attr_id, type_id, flav);
      H5Aclose(attr_id);
      H5Tclose(type_id);
      H5Sclose(attrdat_id);
    }
  
    {
      char *proj = strdup(proj_to_str(corr_thrp.proj));
      hid_t attrdat_id = H5Screate(H5S_SCALAR);
      hid_t type_id = H5Tcopy(H5T_C_S1);
      H5Tset_size(type_id, strlen(proj));
      hid_t attr_id = H5Acreate1(group_id, "Projector", type_id, attrdat_id, H5P_DEFAULT);
      H5Awrite(attr_id, type_id, proj);
      H5Aclose(attr_id);
      H5Tclose(type_id);
      H5Sclose(attrdat_id);
    }
    /* */
  
    hid_t filespace = H5Screate_simple(ND+1, dims, NULL);
    hid_t dataset_id = H5Dcreate(group_id, "corr_x", H5T_NATIVE_DOUBLE, filespace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    filespace = H5Dget_space(dataset_id);
    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, starts, NULL, ldims, NULL);
    hid_t subspace = H5Screate_simple(ND+1, ldims, NULL);
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, subspace, filespace, plist_id, buf);
    H5Dclose(dataset_id);
    H5Sclose(filespace);
    H5Sclose(subspace);
    H5Pclose(plist_id);
    H5Gclose(group_id);
  }
  H5Pclose(lcpl_id);
  H5Gclose(top_id);
  H5Fclose(file_id);
  free(buf);
  return;
}
