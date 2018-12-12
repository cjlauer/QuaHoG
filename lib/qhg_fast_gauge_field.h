#ifndef _QHG_FAST_GAUGE_FIELD_H
#define _QHG_FAST_GAUGE_FIELD_H 1

#include <qhg_types.h>

qhg_fast_gauge_field qhg_fast_gauge_field_init(qhg_lattice *lat);
void qhg_fast_gauge_field_finalize(qhg_fast_gauge_field gf);
void qhg_fast_gauge_field_copy(qhg_fast_gauge_field y, qhg_fast_gauge_field x);
void qhg_fast_gauge_field_id(qhg_fast_gauge_field);
void qhg_fast_gauge_field_import_dir(qhg_fast_gauge_field y, qhg_gauge_field x, int dir);
void qhg_fast_gauge_field_import_dir_daggered(qhg_fast_gauge_field y, qhg_gauge_field x, int dir);
void qhg_fast_gauge_field_multiply(qhg_fast_gauge_field z, qhg_fast_gauge_field x, qhg_fast_gauge_field y);
void qhg_fast_gauge_field_shift(qhg_fast_gauge_field y, qhg_fast_gauge_field x, int dir);

#endif /* _QHG_FAST_GAUGE_FIELD_H */
