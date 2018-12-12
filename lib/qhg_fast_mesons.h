#ifndef _QHG_FAST_MESONS_H
#define _QHG_FAST_MESONS_H 1

#include <qhg_types.h>

qhg_mesons_open_correlator qhg_fast_mesons(qhg_fast_spinor_field sp_1, qhg_fast_spinor_field sp_2);
void qhg_write_mesons_open_correlator(char fname[], qhg_mesons_open_correlator corr, char group[]);
void qhg_mesons_open_correlator_finalize(qhg_mesons_open_correlator corr);


  
#endif /* _QHG_MESONS_H */
