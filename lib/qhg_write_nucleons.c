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
#include <qhg_nucleon_defs.h>

void
qhg_write_nucleons(char fname[], qhg_correlator corr, char group[])
{
  qhg_lattice *lat = corr.lat;
  int proc_id = lat->comms->proc_id;
  int am_io_proc = proc_id == 0 ? 1 : 0;
  size_t site_size = corr.site_size;
  unsigned long int lvol = lat->lvol;
  int *pd = lat->comms->proc_dims;
  int *pc = lat->comms->proc_coords;
  int np = lat->comms->nprocs;
  int *ld = lat->ldims;
  int *d = lat->dims;

  hsize_t starts[ND+3] = {pc[0]*ld[0], pc[1]*ld[1], pc[2]*ld[2], pc[3]*ld[3], 0, 0, 0};
  hsize_t dims[ND+3] = {d[0], d[1], d[2], d[3], NS, NS, 2};
  hsize_t ldims[ND+3] = {ld[0], ld[1], ld[2], ld[3], NS, NS, 2};
  
  hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
  H5Pclose(fapl_id);

  hid_t lcpl_id = H5Pcreate(H5P_LINK_CREATE);
  H5Pset_create_intermediate_group(lcpl_id, 1);  
  hid_t top_id = H5Gcreate(file_id, group, lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
  
  double *buf = qhg_alloc(sizeof(double)*lvol*NS*NS*2);
  _Complex double *c = corr.C;
  /* 
     Hierarchy is /channel/interpolator/ with channel one of ppm or
     pmm (proton or neutron), and with interpolator one of 1-1, 1-2,
     2-1, 2-2 (chi_1, chi_2 combinations)
  */
  for(int iflav=0; iflav<NFLAV; iflav++) {
    char *group_tag;
    asprintf(&group_tag, "%s", flav_tags[iflav]);
    hid_t group1_id = H5Gcreate(top_id, group_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    free(group_tag);

    for(int ichan=0; ichan<NCHAN; ichan++) {
      for(unsigned long int v=0; v<lvol; v++)
	for(int s0=0; s0<NS; s0++)
	  for(int s1=0; s1<NS; s1++) {
	    buf[0 + 2*(s1 + NS*(s0 + NS*v))] = creal(corr.C[NIDX(v, iflav, ichan, s0, s1)]);
	    buf[1 + 2*(s1 + NS*(s0 + NS*v))] = cimag(corr.C[NIDX(v, iflav, ichan, s0, s1)]);
	  }
	  
      asprintf(&group_tag, "%s", chan_tags[ichan]);
      hid_t group2_id = H5Gcreate(group1_id, group_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);      
      free(group_tag);
      /*
	Attributes (metadata) are: 1) the origin (source position) and 2)
	the index order in the file
      */
      hsize_t n = 4;
      hid_t attrdat_id = H5Screate_simple(1, &n, NULL);
      hid_t attr_id = H5Acreate2(group2_id, "Origin", H5T_NATIVE_INT, attrdat_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attr_id, H5T_NATIVE_INT, corr.origin);
      H5Aclose(attr_id);
      H5Sclose(attrdat_id);

      char order[] = "C-order: [t,x,y,z,spin0,spin1,real/imag]\0";
      attrdat_id = H5Screate(H5S_SCALAR);
      hid_t type_id = H5Tcopy(H5T_C_S1);
      H5Tset_size(type_id, strlen(order));
      attr_id = H5Acreate1(group2_id, "IndexOrder", type_id, attrdat_id, H5P_DEFAULT);
      H5Awrite(attr_id, type_id, &order);

      H5Aclose(attr_id);
      H5Tclose(type_id);
      H5Sclose(attrdat_id);
      /* */

      hid_t filespace = H5Screate_simple(ND+3, dims, NULL); 
      hid_t dataset_id = H5Dcreate(group2_id, "corr_x", H5T_NATIVE_DOUBLE, filespace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      hid_t subspace = H5Screate_simple(ND+3, ldims, NULL);
      filespace = H5Dget_space(dataset_id);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, starts, NULL, ldims, NULL);
      hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);      
      herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, subspace, filespace, plist_id, buf);
      H5Dclose(dataset_id);
      H5Sclose(filespace);
      H5Sclose(subspace);
      H5Pclose(plist_id);
      H5Gclose(group2_id);
    }
    H5Gclose(group1_id);
  }
  H5Pclose(lcpl_id);
  H5Gclose(top_id);
  H5Fclose(file_id);
  free(buf);
  return;
}
