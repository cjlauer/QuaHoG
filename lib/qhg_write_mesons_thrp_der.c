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

void
qhg_write_mesons_thrp_der(char fname[], qhg_der_correlator corr, char group[])
{
  qhg_lattice *lat = corr.lat;
  int proc_id = lat->comms->proc_id;
  int am_io_proc = proc_id == 0 ? 1 : 0;
  size_t site_size = corr.site_size;
  int *pd = lat->comms->proc_dims;
  int *pc = lat->comms->proc_coords;
  int np = lat->comms->nprocs;
  int *ld = lat->ldims;
  int *d = lat->dims;

  unsigned long int offset = 0;
  hsize_t starts[ND+2] = {pc[0]*ld[0], pc[1]*ld[1], pc[2]*ld[2], pc[3]*ld[3], 0, 0};
  hsize_t dims[ND+2] = {d[0], d[1], d[2], d[3], 16, 2};
  hsize_t ldims[ND+2] = {ld[0], ld[1], ld[2], ld[3], 16, 2};

  {
    int t0 = corr.origin[0];
    int Lt = dims[0];
    int lt = ldims[0];
    int t_proc_t0 = t0/lt;		/* t-rank of processes with the source time-slice */

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
  
  /*
    Hierarchy is /insertion/ with insertions listed in the
    qhg_nn_thrp_defs.h header file
  */
  for(int id=0; id<(1<<(2*corr.der_order)); id++) {
    if(corr.C[id] != NULL) {
      int dirs[corr.der_order];
      int tmp = id;
      for(int j=corr.der_order-1; j>=0; j--) {
        dirs[j] = tmp % ND;
        tmp /= ND;
      }
      char *group_tag, *tmp_tag;
      asprintf(&group_tag, "der%d:", corr.der_order);
      char c_dirs[]={'0','x','y','z'};
      for(int j=0; j<corr.der_order; j++) {
        asprintf(&tmp_tag, "%sD%c", group_tag, c_dirs[dirs[j]]);
        free(group_tag);
        group_tag = tmp_tag;
      }
      hid_t group_id = H5Gcreate(top_id, group_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      free(group_tag);

      /*
        Attributes (metadata) are:
        1) the origin (source position) [1-dimensional integer array of 4 elements]
        2) the index order in the file [string]
        3) source-sink separation [single integer]
        4) flavor [string]
      */
      {
        hsize_t n = 4;
        hid_t attrdat_id = H5Screate_simple(1, &n, NULL);
        hid_t attr_id = H5Acreate2(group_id, "Origin", H5T_NATIVE_INT, attrdat_id, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_id, H5T_NATIVE_INT, corr.origin);
        H5Aclose(attr_id);
        H5Sclose(attrdat_id);
      }

      {
        char order[] = "C-order: [t,x,y,z,gammas,real/imag]\0";
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
        H5Awrite(attr_id, H5T_NATIVE_INT, &corr.dt);
        H5Aclose(attr_id);
        H5Sclose(attrdat_id);
      }
  
      {
        char *flav;
        switch(corr.flav) {
        case up:
          flav = strdup("up");
          break;
        case dn:
          flav = strdup("dn");
          break;
	case strange:
	  flav = strdup("strange");
	  break;
	default:
	  fprintf(stderr, " Flavor %s is not supported\n", corr.flav);
	  exit(3);
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
  
      int lt = lat->ldims[0];
      unsigned long int lv3 = lat->lvol/(unsigned long int)lt;
      double *buf = ((double*) corr.C[id]) + offset*lv3*16*2;
      hid_t filespace = H5Screate_simple(ND+2, dims, NULL);
      hid_t dataset_id = H5Dcreate(group_id, "corr_x", H5T_NATIVE_DOUBLE, filespace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      filespace = H5Dget_space(dataset_id);
      hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, starts, NULL, ldims, NULL);
      hid_t subspace = H5Screate_simple(ND+2, ldims, NULL);
      herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, subspace, filespace, plist_id, buf);
      H5Dclose(dataset_id);
      H5Sclose(filespace);
      H5Sclose(subspace);
      H5Pclose(plist_id);
      H5Gclose(group_id);
    }
  }
  H5Pclose(lcpl_id);
  H5Gclose(top_id);
  H5Fclose(file_id);

  return;
}

void
qhg_write_mesons_averaged_thrp_der(char fname[], qhg_der_correlator corr, char group[])
{
  qhg_lattice *lat = corr.lat;
  int proc_id = lat->comms->proc_id;
  int am_io_proc = proc_id == 0 ? 1 : 0;
  size_t site_size = corr.site_size;
  int *pd = lat->comms->proc_dims;
  int *pc = lat->comms->proc_coords;
  int np = lat->comms->nprocs;
  int *ld = lat->ldims;
  int *d = lat->dims;
  int ncorr = corr.ncorr;
  int der_order = corr.der_order;

  unsigned long int offset = 0;
  hsize_t starts[ND+1] = {pc[0]*ld[0], pc[1]*ld[1], pc[2]*ld[2], pc[3]*ld[3], 0};
  hsize_t dims[ND+1] = {d[0], d[1], d[2], d[3], 2};
  hsize_t ldims[ND+1] = {ld[0], ld[1], ld[2], ld[3], 2};

  {
    int t0 = corr.origin[0];
    int Lt = dims[0];
    int lt = ldims[0];
    int t_proc_t0 = t0/lt;		/* t-rank of processes with the source time-slice */

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
  
  char *c_dirs;
  if(der_order == 2){

    if(ncorr == 3){
      char c_dirs_tmp[3*3] = { '0','x','y',
			       '0','x','z',
			       '0','y','z' };
      c_dirs = c_dirs_tmp;

    } 
    else {
      printf("ERROR: Number of averaged derivative combination %d not supported for derivative order %d.\n", ncorr, der_order);
      return;
    }

  }
  else if(der_order == 3){

    if(ncorr == 6){
      char c_dirs_tmp[6*4] = { '0','0','x','x',
			       '0','0','y','y',
			       '0','0','z','z',
			       'x','x','y','y',
			       'x','x','z','z',
			       'y','y','z','z' };
      c_dirs = c_dirs_tmp;
    } 
    else if(ncorr == 7){
      char c_dirs_tmp[7*4] = { '0','0','x','x',
			       '0','0','y','y',
			       '0','0','z','z',
			       'x','x','y','y',
			       'x','x','z','z',
			       'y','y','z','z',
			       '0','x','y','z' };
      c_dirs = c_dirs_tmp;
    } 
    else {
      printf("ERROR: Number of averaged derivative combination %d not supported for derivative order %d.\n", ncorr, der_order);
      return;
    }

  } 
  else {

    printf("ERROR: derivative order %d not supported by write function\n", der_order);
    return;

  }

  for( int c=0; c<ncorr; c++) {
    if(corr.C[c] != NULL) {

      char *group_tag, *tmp_tag;
      asprintf(&group_tag, "der%d:g%c", der_order, c_dirs[(der_order+1)*c]);
      for(int j=1; j<=der_order; j++) {
	asprintf(&tmp_tag, "%sD%c", group_tag, c_dirs[(der_order+1)*c+j]);
	free(group_tag);
	group_tag = tmp_tag;
      }
      hid_t group_id = H5Gcreate(top_id, group_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      free(group_tag);

      /*
	Attributes (metadata) are:
	1) the origin (source position) [1-dimensional integer array of 4 elements]
	2) the index order in the file [string]
	3) source-sink separation [single integer]
	4) flavor [string]
      */

      {
        hsize_t n = 4;
        hid_t attrdat_id = H5Screate_simple(1, &n, NULL);
        hid_t attr_id = H5Acreate2(group_id, "Origin", H5T_NATIVE_INT, attrdat_id, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_id, H5T_NATIVE_INT, corr.origin);
        H5Aclose(attr_id);
        H5Sclose(attrdat_id);
      }

      {
        char order[] = "C-order: [t,x,y,z,gammas,real/imag]\0";
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
        H5Awrite(attr_id, H5T_NATIVE_INT, &corr.dt);
        H5Aclose(attr_id);
        H5Sclose(attrdat_id);
      }
  
      {
        char *flav;
        switch(corr.flav) {
        case up:
          flav = strdup("up");
          break;
        case dn:
          flav = strdup("dn");
          break;
	case strange:
	  flav = strdup("strange");
	  break;
	default:
	  fprintf(stderr, "Flavor %s is not supported\n", corr.flav);
	  exit(3);
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
  
      int lt = lat->ldims[0];
      unsigned long int lv3 = lat->lvol/(unsigned long int)lt;
      double *buf = ((double*) corr.C[c]) + offset*lv3*2;
      hid_t filespace = H5Screate_simple(ND+1, dims, NULL);
      hid_t dataset_id = H5Dcreate(group_id, "corr_x", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
  }
  
  H5Pclose(lcpl_id);
  H5Gclose(top_id);
  H5Fclose(file_id);

  return;
}
