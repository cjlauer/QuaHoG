#ifndef _QHG_FAST_BARYONS_H
#define _QHG_FAST_BARYONS_H 1

#include <qhg_types.h>

qhg_baryons_open_correlator qhg_fast_baryons(qhg_fast_spinor_field sp_1, qhg_fast_spinor_field sp_2, qhg_fast_spinor_field sp_3);
void qhg_write_baryons_open_correlator(char fname[], qhg_baryons_open_correlator corr, char group[]);
void qhg_baryons_open_correlator_finalize(qhg_baryons_open_correlator corr);


  
#endif /* _QHG_BARYONS_H */
