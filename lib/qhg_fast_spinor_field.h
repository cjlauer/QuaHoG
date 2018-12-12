#ifndef _QHG_FAST_SPINOR_FIELD_H
#define _QHG_FAST_SPINOR_FIELD_H 1

#include <qhg_types.h>

qhg_fast_spinor_field qhg_fast_spinor_field_init(qhg_lattice *, enum qhg_fermion_bc_time);
void qhg_fast_spinor_field_finalize(qhg_fast_spinor_field );
void qhg_fast_spinor_field_copy(qhg_fast_spinor_field, qhg_fast_spinor_field);
void qhg_fast_spinor_field_import(qhg_fast_spinor_field, qhg_spinor_field[NS*NC]);
void qhg_fast_spinor_field_import_propagator(qhg_fast_spinor_field y, qhg_spinor_field x, int cs0);
void qhg_fast_spinor_field_remove_bc_lower(qhg_fast_spinor_field x, int tsrc);
void qhg_fast_spinor_field_remove_bc_higher(qhg_fast_spinor_field x, int tsrc);
double qhg_fast_spinor_field_normsq(qhg_fast_spinor_field);
void qhg_fast_spinor_field_axpy(double,qhg_fast_spinor_field,qhg_fast_spinor_field);
void qhg_fast_spinor_field_swap(qhg_fast_spinor_field*, qhg_fast_spinor_field*);
void qhg_fast_spinor_field_shift(qhg_fast_spinor_field,qhg_fast_spinor_field, int);
void qhg_fast_spinor_field_multiply_G_G(qhg_fast_spinor_field, qhg_fast_spinor_field, qhg_fast_spinor_field);
void qhg_fast_spinor_field_trace_multiply_U_G(qhg_fast_gauge_field U, qhg_fast_spinor_field y,
                                              qhg_correlator corr, int id, double factor);
void qhg_fast_spinor_field_trace_multiply_Udag_G(qhg_fast_gauge_field U, qhg_fast_spinor_field y,
                                                 qhg_correlator corr, int id, double factor);
void qhg_fast_spinor_field_trace_multiply_G_U_G(qhg_fast_spinor_field x, qhg_fast_gauge_field U, qhg_fast_spinor_field y,
                                              qhg_correlator corr, int id, double factor);
void qhg_fast_spinor_field_trace_multiply_G_Udag_G(qhg_fast_spinor_field x, qhg_fast_gauge_field U, qhg_fast_spinor_field y,
                                              qhg_correlator corr, int id, double factor);
#endif /* _QHG_FAST_SPINOR_FIELD_H */
