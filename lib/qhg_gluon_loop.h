#ifndef _QHG_GLUON_LOOP_H
#define _QHG_GAUGE_FIELD_H 1

#include <qhg_types.h>

qhg_gluon_loop qhg_gluon_loop_init(qhg_lattice *);
void qhg_gluon_loop_copy(qhg_gluon_loop, qhg_gluon_loop);
static double gluon_plaq(_Complex double *, unsigned long int **, unsigned long int, int, int);
qhg_gluon_loop qhg_calculate_gluon_loop(qhg_gluon_loop, qhg_gauge_field);
void qhg_gluon_loop_finalize(qhg_gluon_loop);
  
#endif /* _QHG_GLUON_LOOP_H */
