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

qhg_mesons_open_correlator
qhg_fast_mesons(qhg_fast_spinor_field sp_1, qhg_fast_spinor_field sp_2)
{
  qhg_lattice *lat = sp_1.lat; 
  int lt = lat->ldims[0];
  qhg_mesons_open_correlator corr;
  corr.C = qhg_alloc(lt*NS*NS*NS*NS*2*sizeof(double));
  memset(corr.C, '\0', lt*NS*NS*NS*NS*2*sizeof(double));
  corr.lat = lat; 
  unsigned long int lv3 = lat->lv3;

#ifdef QHG_OMP
#pragma omp parallel for
#endif
  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++)
        for(int s2=0; s2<NS; s2++)
          for(int s3=0; s3<NS; s3++)
            for(int c0=0; c0<NC; c0++)
              for(int c1=0; c1<NC; c1++)
                for(unsigned long int v=0; v<lv3; v++) {
                  corr.C[t][s0][s1][s2][s3][0] += sp_1.field[t][s0][s1][c0][c1][0][v] * sp_2.field[t][s2][s3][c0][c1][0][v] + sp_1.field[t][s0][s1][c0][c1][1][v] * sp_2.field[t][s2][s3][c0][c1][1][v];
                  corr.C[t][s0][s1][s2][s3][1] += sp_1.field[t][s0][s1][c0][c1][0][v] * sp_2.field[t][s2][s3][c0][c1][1][v] - sp_1.field[t][s0][s1][c0][c1][1][v] * sp_2.field[t][s2][s3][c0][c1][0][v];
                }


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
    MPI_Reduce( MPI_IN_PLACE, corr.C, lt*NS*NS*NS*NS*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    corr.write = 1;
  } else {
    MPI_Reduce( corr.C, NULL, lt*NS*NS*NS*NS*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    corr.write = 0;
  }

  //MPI_Comm_free(&tcomm);
  return corr;
}

void
qhg_mesons_open_correlator_finalize(qhg_mesons_open_correlator corr)
{
  free(corr.C);
  corr.lat = NULL;
  corr.write = 0;
}
  
void
qhg_write_mesons_open_correlator(char fname[], qhg_mesons_open_correlator corr, char group[])
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

  int ndims = 6;
  hsize_t starts[6] = {pc[0]*ld[0], 0, 0, 0, 0, 0};
  hsize_t dims[6] = {d[0], 4, 4, 4, 4, 2};
  hsize_t ldims[6] = {ld[0], 4, 4, 4, 4, 2};
  
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

  char order[] = "C-order: [t,s0,s1,s3,s4,re-re/re-im/im-re/im-im]\0";
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
