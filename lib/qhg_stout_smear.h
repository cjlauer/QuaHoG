#ifndef _QHG_STOUT_SMEAR
#define _QHG_STOUT_SMEAR 1

#include <qhg_types.h>

double qhg_xi0(double);
void qhg_exp(_Complex double *);
void qhg_stout_smear_iter(qhg_gauge_field, qhg_gauge_field, double);
void qhg_stout_smear(qhg_gauge_field, qhg_gauge_field, double, int);

#endif /* _QHG_STOUT_SMEAR */
