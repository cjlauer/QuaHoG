#ifndef _QHG_XCHANGE_SPINOR_H
#define _QHG_XCHANGE_SPINOR_H 1

#include <qhg_types.h>

void get_boundary_spinor(_Complex double *bnd, int dir, _Complex double *field, qhg_lattice *lat);
void qhg_xchange_spinor(qhg_spinor_field);

#endif /* _QHG_XCHANGE_SPINOR_H */
