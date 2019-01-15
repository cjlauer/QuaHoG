#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

#include <complex.h>
#include <qhg_idx.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_gluon_loop.h>
#include <qhg_prop_gammas.h>
#include <qhg_prop_ops.h>
#include <qhg_io_utils.h>
#include <qhg_meson_defs.h>

void
qhg_write_gluon_loops(char fname[], qhg_gluon_loop gl, char group[])
{
  qhg_lattice *lat = gl.lat;
  int proc_id = lat->comms->proc_id;
  int am_io_proc = proc_id == 0 ? 1 : 0;
  unsigned long int lvol = lat->lvol;
  int *pd = lat->comms->proc_dims;
  int *pc = lat->comms->proc_coords;
  int np = lat->comms->nprocs;
  int *ld = lat->ldims;
  int *d = lat->dims;

  hsize_t starts[ND+1] = {pc[0]*ld[0], pc[1]*ld[1], pc[2]*ld[2], pc[3]*ld[3], 0};
  hsize_t dims[ND+1] = {d[0], d[1], d[2], d[3], 2};
  hsize_t ldims[ND+1] = {ld[0], ld[1], ld[2], ld[3], 2};
  
  hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
  H5Pclose(fapl_id);

  hid_t lcpl_id = H5Pcreate(H5P_LINK_CREATE);
  H5Pset_create_intermediate_group(lcpl_id, 1);  
  hid_t top_id = H5Gcreate(file_id, group, lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
  
  double *buf = qhg_alloc(sizeof(double)*lvol*2);
  _Complex double *loop = gl.loop;
  for(unsigned long int v=0; v<lvol; v++) {
    buf[0 + 2*v] = creal(loop[v]);
    buf[1 + 2*v] = cimag(loop[v]);
  }

  /*
    Attributes (metadata) are: 
    1) the index order in the file
  */
  char order[] = "C-order: [t,x,y,z,real/imag]\0";
  hid_t attrdat_id = H5Screate(H5S_SCALAR);
  hid_t type_id = H5Tcopy(H5T_C_S1);
  H5Tset_size(type_id, strlen(order));
  hid_t attr_id = H5Acreate1(top_id, "IndexOrder", type_id, attrdat_id, H5P_DEFAULT);
  H5Awrite(attr_id, type_id, &order);

  H5Aclose(attr_id);
  H5Tclose(type_id);
  H5Sclose(attrdat_id);
  /* */
      
  hid_t filespace = H5Screate_simple(ND+1, dims, NULL); 
  hid_t dataset_id = H5Dcreate(top_id, "gluon_loop_x", H5T_NATIVE_DOUBLE, filespace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t subspace = H5Screate_simple(ND+1, ldims, NULL);
  filespace = H5Dget_space(dataset_id);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, starts, NULL, ldims, NULL);
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);      
  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, subspace, filespace, plist_id, buf);
  H5Dclose(dataset_id);
  H5Sclose(filespace);
  H5Sclose(subspace);
  H5Pclose(plist_id);

  H5Pclose(lcpl_id);
  H5Gclose(top_id);
  H5Fclose(file_id);
  return;
}
