#ifndef _QHG_NN_THRP_DER_H
#define _QHG_NN_THRP_DER_H 1

#include <qhg_types.h>
#include <stdbool.h>

qhg_der_correlator
qhg_nn_thrp_der(qhg_spinor_field fwd[NS*NC], qhg_spinor_field bwd[NS*NC], qhg_gauge_field gf,
                int source_coords[ND], qhg_thrp_nn_sink_params thrp_sink, int der_order, bool (*to_skip)(int*),
                bool rev_trick, bool perm_trick);

#endif /* _QHG_NN_THRP_H */
