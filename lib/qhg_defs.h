#ifndef _QHG_DEFS_H
#define _QHG_DEFS_H 1
/*
  Number of colors
  Number of spins
  Number of dimensions
*/
#define NC 3
#define NS 4
#define ND 4
#define SPINOR_EDGES

#define QHG_MEMALIGN 64
#define QHG_FAST_TYPE float
typedef QHG_FAST_TYPE afloat __attribute__ ((aligned (QHG_MEMALIGN)));

#define NEPS 6
extern int QHG_EPS[NEPS][NC+1];

#define MAX(a,b) (a) > (b) ? (a) : (b)
#define MIN(a,b) (a) < (b) ? (a) : (b)

#endif /* _QHG_DEFS_H */
