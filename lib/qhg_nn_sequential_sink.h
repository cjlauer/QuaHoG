#ifndef _QHG_NN_SEQUENTIAL_SINKS_H
#define _QHG_NN_SEQUENTIAL_SINKS_H 1

#include <qhg_types.h>

void qhg_nn_sequential_sink_d(qhg_spinor_field [], qhg_spinor_field [], int, qhg_thrp_nn_sink_params);
void qhg_nn_sequential_sink_u(qhg_spinor_field [], qhg_spinor_field [], qhg_spinor_field [], int, qhg_thrp_nn_sink_params);
void qhg_mesons_sequential_sink(qhg_spinor_field [], qhg_spinor_field [], int, qhg_thrp_nn_sink_params);
void qhg_phase_sequential_sink(qhg_spinor_field [], qhg_spinor_field [], int mom[3], int source_coords[4], int);

  
#endif /* _QHG_NN_SEQUENTIAL_SINKS_H */
