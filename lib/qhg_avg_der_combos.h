#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include <complex.h>
#include <qhg_idx.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_correlator.h>
#include <qhg_correlator_shift.h>
#include <qhg_spinor_field.h>
#include <qhg_fast_spinor_field.h>
#include <qhg_fast_gauge_field.h>
#include <qhg_xchange_spinor.h>
#include <qhg_xchange_gauge.h>
#include <qhg_prop_gammas.h>
#include <qhg_prop_ops.h>
#include <qhg_su3_ops.h>
#include <qhg_io_utils.h>
#include <qhg_nn_thrp_der.h>
#include <math.h>
#define VG(v, g) (v*NS*NS+g)
#define DD(mu,nu) (ND*mu+nu)
#define DDD(mu,nu,rho) (ND*ND*mu+ND*nu+rho)

/* 
   Calculates average of 2nd derivative correlators of the form
   O_ijk = g_i D_j D_k over permutations of two indices and 0 
*/

void qhg_avg_der_combos_O_0_i_j(_Complex double *, _Complex double **, int, int, size_t);

/*
  Calculate average of 3rd derivative correlators of the form 
  O_{mu mu nu nu} = g_mu D_mu D_nu D_nu over permutations of two indices
*/
void qhg_avg_der_combos_O_mu_mu_nu_nu(_Complex double *, _Complex double **, int, int, size_t);

/*
  Calculate average of 3rd derivative correlators of the form 
  O_{mu nu rho sig} g_mu D_nu D_rho D_sig over permutations 
  of four indices without any pairs of equal indices
*/
void qhg_avg_der_combos_O_mu_nu_rho_sig(_Complex double *, _Complex double **, int, int, int, int, size_t);

qhg_der_correlator qhg_avg_der_combos_der2(qhg_der_correlator);

qhg_der_correlator qhg_avg_der_combos_der3(qhg_der_correlator, int, int, int);

qhg_der_correlator qhg_avg_der_combos(qhg_der_correlator, int, int, int);
