#ifndef _QHG_IDX_H
#define _QHG_IDX_H 1
#include <qhg_defs.h>
/*
 * Helper macros for indexing
 */

/*
 * Single color-spin index from [spin,color] index
 */
#define CS(s,c) ((c)+(s)*NC)
#define CS2S(cs) ((int) cs/NC)
#define CS2C(cs) ((int) cs%NC)

/*
 * Single index from [color0,color1] index
 */
#define CC(i,j) ((j)+(i)*NC)

/*
 * Returns {t, x, y, z} coordinate corresponding to single index v,
 * when dimensions are {d[0], d[1], d[2], d[3]}.  To be used in
 * initializer, e.g. int x[ND] = CO(v, dims);
 */
#define CO(v,d) {((v) / (d[1]*d[2]*d[3])),				\
		 ((v) / (d[3]*d[2])) % d[1],				\
		 ((v) / d[3]) % d[2],					\
		 ((v) % d[3])						};

/*
 * Returns single index of coordinates {co[0], co[1], co[2], co[3]}
 * when dimensions are {d[0], d[1], d[2], d[3]}
 */
#define IDX(co,d) ((unsigned long int)co[3] +				\
		   (unsigned long int)d[3]*((unsigned long int)co[2] +	\
					    (unsigned long int)d[2]*((unsigned long int)co[1] + \
								     (unsigned long int)d[1]*(unsigned long int)co[0])))

/*
 * Same as above but for three dimensions
 * 
 */
#define IDX3(co,d) ((unsigned long int)co[2] + \
		    (unsigned long int)d[2]*((unsigned long int)co[1] + \
					     (unsigned long int)d[1]*(unsigned long int)co[0]))

#endif /* _QHG_IDX_H */
