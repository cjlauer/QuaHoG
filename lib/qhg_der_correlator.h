#include <string.h>

#include <qhg_alloc.h>
#include <qhg_defs.h>
#include <qhg_types.h>

qhg_der_correlator
qhg_der_correlator_init(size_t, qhg_lattice *, qhg_mom_list *, int);

qhg_der_correlator
qhg_averaged_der_correlator_init(size_t, qhg_lattice *, qhg_mom_list *, int);

qhg_der_correlator
qhg_der_correlator_copy(qhg_der_correlator);

void
qhg_correlator_finalize(qhg_der_correlator);
