#ifndef _QHG_NN_THRP_DEFS_DER2_H
#define _QHG_NN_THRP_DEFS_DER2_H 1
#define NCHAN 16	// number of gamma channels
#define SITE_SIZE (NCHAN*ND*ND)
#define CHAN_IDX(v, ch) (ch + v*NCHAN)
#define TOT_IDX(v,mu,nu,ch) (ch + nu*NCHAN + mu*ND*NCHAN + v*NCHAN*ND*ND)

static char chan_tags[NCHAN][256] = {
  "=der2:one", "=der2:g5",
  "=der2:g0", "=der2:gx", "=der2:gy", "=der2:gz", 
  "=der2:g5g0", "=der2:g5gx", "=der2:g5gy", "=der2:g5gz", 
  "=der2:g5si0x", "=der2:g5si0y", "=der2:g5si0z", 
  "=der2:g5sixy", "=der2:g5sixz", "=der2:g5siyz"
};

static char der_tags[ND][ND][256] = {
  { "D0D0=", "D0Dx=", "D0Dy=", "D0Dz=", },
  { "DxD0=", "DxDx=", "DxDy=", "DxDz=", },
  { "DyD0=", "DyDx=", "DyDy=", "DyDz=", },
  { "DzD0=", "DzDx=", "DzDy=", "DzDz=", },
};

static enum {
  one, g5, g0, gx, gy, gz,
  g5g0, g5gx, g5gy, g5gz,
  g5si0x, g5si0y, g5si0z, g5sixy, g5sixz, g5siyz,
} channels;

#endif /* _QHG_NN_THRP_DEFS_H */
