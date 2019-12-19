#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <mpi.h>
#include <mxml.h>
#include <lime.h>

#include <qhg_types.h>
#include <qhg_idx.h>
#include <qhg_defs.h>
#include <qhg_alloc.h>
#include <qhg_io_utils.h>

static struct {
  int dims[ND];
  int precision;
  n_uint64_t binary_data_offset;
} conf_info;

static int
get_int_from_xml(char s[], char *buf)
{
  mxml_node_t *tree = mxmlLoadString(NULL, buf, MXML_TEXT_CALLBACK);  
  char *path;
  asprintf(&path, "ildgFormat/%s", s);
  mxml_node_t *node = mxmlFindPath(tree, path);
  const char *p = mxmlGetText(node, NULL);
  int r = atoi(p);
  free(path);
  mxmlDelete(node);    
  mxmlDelete(tree);      
  return r;
}

static void
handle_ildg_format(LimeReader *reader)
{
  n_uint64_t nbytes = limeReaderBytes(reader);
  char *buf = malloc((size_t)nbytes+1);
  n_uint64_t bytes_read = nbytes;
  int status = limeReaderReadData(buf, &bytes_read, reader);
  conf_info.precision = get_int_from_xml("precision", buf);
  conf_info.dims[1] = get_int_from_xml("lx", buf);
  conf_info.dims[2] = get_int_from_xml("ly", buf);
  conf_info.dims[3] = get_int_from_xml("lz", buf);
  conf_info.dims[0] = get_int_from_xml("lt", buf);
  free(buf);  
  return;
}

static void
handle_ildg_binary_data(LimeReader *reader)
{
  n_uint64_t bytes_in_record = limeReaderBytes(reader);
  n_uint64_t expected_bytes =
    (n_uint64_t)(conf_info.precision/8) *
    (n_uint64_t)conf_info.dims[0] *
    (n_uint64_t)conf_info.dims[1] *
    (n_uint64_t)conf_info.dims[2] *
    (n_uint64_t)conf_info.dims[3] *
    (n_uint64_t)(NC * NC * ND * 2);
  if(expected_bytes != bytes_in_record) {
    fprintf(stderr, "expected %llu bytes in ildg-binary-data, got %llu\n",
	    expected_bytes, bytes_in_record);
    exit(2);
  }
  return;
}

static void
get_file_types(MPI_Datatype *etype, MPI_Datatype *ftype, MPI_Datatype dtype, qhg_lattice *lat)
{
  MPI_Type_contiguous(2*NC*NC*ND, dtype, etype);
  MPI_Type_commit(etype);
  
  int ldims[] = {
    lat->ldims[0], lat->ldims[3], lat->ldims[2], lat->ldims[1]
  };
  int dims[] = {
    lat->dims[0], lat->dims[3], lat->dims[2], lat->dims[1]
  };

  int *proc_coords = lat->comms->proc_coords;
  int starts[ND];
  starts[0] = proc_coords[0]*ldims[0];
  starts[3] = proc_coords[1]*ldims[3];
  starts[2] = proc_coords[2]*ldims[2];
  starts[1] = proc_coords[3]*ldims[1];
  int ierr = MPI_Type_create_subarray(ND, dims, ldims, starts, MPI_ORDER_C, *etype, ftype);
  MPI_Type_commit(ftype);
  return;
}
 
void
qhg_read_gauge_field_ildg(qhg_gauge_field g, char fname[])
{
  qhg_lattice *lat = g.lat;
  int *dims = lat->dims;
  int am_io_proc = lat->comms->proc_id == 0 ? 1 : 0;
  FILE *fp = qhg_fopen(fname, "r");
  LimeReader *reader = limeCreateReader(fp);
  if(reader == NULL) {
    if(am_io_proc) {
      fprintf(stderr, " %s: error opening as lime file\n", fname);      
    }
    exit(2);
  }

  int status = limeReaderNextRecord(reader);
  while(status != LIME_EOF) {
    char *type = limeReaderType(reader);

    if(strcmp(type, "ildg-format") == 0) {
      handle_ildg_format(reader);
    }
    if(strcmp(type, "ildg-binary-data") == 0) {
      handle_ildg_binary_data(reader);
      conf_info.binary_data_offset = ftell(fp);
      break;
    }
    status = limeReaderNextRecord(reader);    
  }

  if((conf_info.dims[0] != dims[0]) ||
     (conf_info.dims[3] != dims[3]) ||
     (conf_info.dims[2] != dims[2]) ||
     (conf_info.dims[1] != dims[1])) {
    if(am_io_proc) {
      fprintf(stderr, " %s: dimensions do not match\n", fname);
    }
    exit(2);
  }

  if(conf_info.precision != 32 &
     conf_info.precision != 64) {
    if(am_io_proc) {
      fprintf(stderr, " %s: precision neither 32 nor 64\n", fname);      
    }
    exit(2);
  }
  
  limeDestroyReader(reader);
  fclose(fp);

  size_t elem_bytes = conf_info.precision == 32 ? 4 : 8;
  MPI_Datatype dtype = elem_bytes == 4 ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Datatype elemtype, filetype;
  get_file_types(&elemtype, &filetype, dtype, lat);
  
  MPI_File fhandle;
  int err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhandle);
  if(err != 0) {
    if(am_io_proc) {
      fprintf(stderr, "error opening file in %s\n", __func__);      
    }
    MPI_Abort(MPI_COMM_WORLD, -2);
  }
  MPI_File_set_view(fhandle, (MPI_Offset)conf_info.binary_data_offset, elemtype,
		    filetype, "native", MPI_INFO_NULL);

  unsigned long int lvol = lat->lvol;
  void *buffer = qhg_alloc(lvol*2*elem_bytes*NC*NC*ND);
  MPI_File_read_all(fhandle, buffer, (int)lvol, elemtype, MPI_STATUS_IGNORE);
  if(!qhg_is_bigendian()) {
    if(elem_bytes == 4) {
      qhg_byte_swap_float(buffer, lvol*NC*NC*ND*2);
    } else {
      qhg_byte_swap_double(buffer, lvol*NC*NC*ND*2);
    }
  }
  
  for(unsigned long int v0=0; v0<lvol; v0++) {
    int *ldims0 = lat->ldims;    
    int co0[] = CO(v0, ldims0);
    int co1[] = {
      co0[0], co0[3], co0[2], co0[1]
    };
    int ldims1[] = {
      ldims0[0], ldims0[3], ldims0[2], ldims0[1]
    };    
    unsigned long int v1 = IDX(co1, ldims1);
    int mu[] = {3, 0, 1, 2};
    for(int d=0; d<ND; d++) 
      for(int c0=0; c0<NC; c0++)
	for(int c1=0; c1<NC; c1++) {
	  unsigned long int x0 = (unsigned long int)c1 +
	    NC*((unsigned long int)c0 +
		NC*((unsigned long int)d + (unsigned long int)v0*ND));
	  unsigned long int x1 = (unsigned long int)c1 +
	    NC*((unsigned long int)c0 + NC*((unsigned long int)mu[d] + (unsigned long int)v1*ND));
	  double re, im;
	  if(elem_bytes == 4) {
	    re = ((float *) buffer)[x1*2 + 0];
	    im = ((float *) buffer)[x1*2 + 1];
	  }
	  if(elem_bytes == 8) {
	    re = ((double *) buffer)[x1*2 + 0];
	    im = ((double *) buffer)[x1*2 + 1];
	  }
	  g.field[x0] = re + _Complex_I*im;
	}
  }
  free(buffer);    

  MPI_File_close(&fhandle);
  MPI_Type_free(&elemtype);
  MPI_Type_free(&filetype);
  
  return;
}
