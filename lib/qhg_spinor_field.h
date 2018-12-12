#ifndef _QHG_SPINOR_FIELD_H
#define _QHG_SPINOR_FIELD_H 1

#include <qhg_types.h>

qhg_spinor_field qhg_spinor_field_init(qhg_lattice *, enum qhg_fermion_bc_time);
void qhg_spinor_field_finalize(qhg_spinor_field );
void qhg_spinor_field_copy(qhg_spinor_field, qhg_spinor_field);
double qhg_spinor_field_normsq(qhg_spinor_field);
void qhg_spinor_field_shift(qhg_spinor_field prop, int shifts_in[ND]);
void qhg_spinor_field_axpy(double, qhg_spinor_field, qhg_spinor_field);

#endif /* _QHG_SPINOR_FIELD_H */
