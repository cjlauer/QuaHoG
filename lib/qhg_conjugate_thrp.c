#include <complex.h>
#include <qhg_idx.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_correlator.h>
#include <qhg_xchange_spinor.h>
#include <qhg_xchange_gauge.h>
#include <qhg_prop_gammas.h>
#include <qhg_prop_ops.h>
#include <qhg_su3_ops.h>
#include <qhg_io_utils.h>
#include <qhg_nn_thrp_defs.h>
#include <qhg_nn_thrp_defs.h>
#include <qhg_conjugate_thrp.h>

// Function to take the complex conjugate of the three-point data
void qhg_conjugate_thrp(qhg_correlator corr, int sign) {

  qhg_lattice *lat = corr.lat;
  unsigned long int lvol = lat->lvol;
  
  // The sign depends on the conjugate of the insertion current's gamma structure

  // Local
  int sign_loc[NLOC] = { 
    // Scalar
    1,
    // Pseudo-scalar
    1,
    // Vector
    -1,
    -1,
    -1,
    -1,
    // Axial
    1,
    1,
    1,
    1,
    // Tensor
    -1,
    -1,
    -1,
    -1,
    -1,
    -1,
  };

  // Noether
  int sign_noe = -1;

  // Vector derivative
  int sign_vder =  1;

  // Axial derivative, symmetric
  int sign_adsy = -1;

  // Axial derivative, anti-symmetric
  int sign_adas = -1;
    // Tensor derivative, symmetric
  int sign_tder = 1;

#ifdef QHG_OMP
#pragma omp parallel for
#endif
  for(unsigned long int v=0; v<lvol; v++) {

    // Local
    for(int i=0; i<NLOC; i++) {
      corr.C[TIDX(v, i)] = 1.0 * sign * sign_loc[i] * conj( corr.C[TIDX(v,i)] );
    }
    // Noether
    for(int i=NLOC; i<NLOC+NNOE; i++) {
      corr.C[TIDX(v, i)] = 1.0 * sign * sign_noe * conj( corr.C[TIDX(v,i)] );
    }
    // Vector derivative
    for(int i=NLOC+NNOE; i<NLOC+NNOE+NVDER; i++) {
      corr.C[TIDX(v, i)] = 1.0 * sign * sign_vder * conj( corr.C[TIDX(v,i)] );
    }      
    /* Axial derivative, symmetric */
    for(int i=NLOC+NNOE+NVDER; i<NLOC+NNOE+NVDER+NADSY; i++) {
      corr.C[TIDX(v, i)] = 1.0 * sign * sign_adsy * conj( corr.C[TIDX(v,i)] );
    }
    /* Axial derivative, anti-symmetric */
    for(int i=NLOC+NNOE+NVDER+NADSY; i<NLOC+NNOE+NVDER+NADSY+NADAS; i++) {
      corr.C[TIDX(v, i)] = 1.0 * sign * sign_adas * conj( corr.C[TIDX(v,i)] );
    }
    /* Tensor derivative, symmetric */
    for(int i=NLOC+NNOE+NVDER+NADSY+NADAS; i<NLOC+NNOE+NVDER+NADSY+NADAS+NTDER; i++) {
      corr.C[TIDX(v, i)] = 1.0 * sign * sign_tder * conj( corr.C[TIDX(v,i)] );
    }
  } // -v
}

