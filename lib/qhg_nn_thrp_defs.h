#ifndef _QHG_NN_THRP_DEFS_H
#define _QHG_NN_THRP_DEFS_H 1
#define NLOC 16			// number of local operator directions
#define NNOE 4			// number of noether operator directions
#define NVDER ((4*(4+1))/2)	// number of vector-derivative directions
#define NADSY ((4*(4+1))/2)	// number of axial-derivative directions (symmetric)
#define NADAS ((4*(4-1))/2)	// number of axial-derivative directions (anti-symmetric)
#define NTDER (4*((4+1)*4)/2-4) // number of tensor-derivative directions [ (si_{mu,nu} D_ku + si_{mu,ku} D_nu) with nu <= ku and ommit ku == mu == nu which is zero]
#define NCHAN (NLOC+NNOE+NVDER+NADSY+NADAS+NTDER)
#define SITE_SIZE (NCHAN)
#define TIDX(v, ch) (ch + v*NCHAN)

static char chan_tags[NCHAN][256] = {
  "=loc:1=\0", "=loc:g5=\0",
  "=loc:g0=\0", "=loc:gx=\0", "=loc:gy=\0", "=loc:gz=\0",
  "=loc:g5g0=\0", "=loc:g5gx=\0", "=loc:g5gy=\0", "=loc:g5gz=\0",
  "=loc:g5si0x=\0", "=loc:g5si0y=\0", "=loc:g5si0z=\0", "=loc:g5sixy=\0", "=loc:g5sixz=\0", "=loc:g5siyz=\0",
  //
  "=noe:g0=\0", "=noe:gx=\0", "=noe:gy=\0", "=noe:gz=\0",
  //
  "=der:g0D0:sym=\0", "=der:gxD0:sym=\0", "=der:gyD0:sym=\0", "=der:gzD0:sym=\0",
  "=der:gxDx:sym=\0", "=der:gyDx:sym=\0", "=der:gzDx:sym=\0",
  "=der:gyDy:sym=\0", "=der:gzDy:sym=\0",
  "=der:gzDz:sym=\0",
  //
  "=der:g5g0D0:sym=\0", "=der:g5gxD0:sym=\0", "=der:g5gyD0:sym=\0", "=der:g5gzD0:sym=\0",
  "=der:g5gxDx:sym=\0", "=der:g5gyDx:sym=\0", "=der:g5gzDx:sym=\0",
  "=der:g5gyDy:sym=\0", "=der:g5gzDy:sym=\0",
  "=der:g5gzDz:sym=\0",
  //
  "=der:g5gxD0:asy=\0", "=der:g5gyD0:asy=\0", "=der:g5gzD0:asy=\0",
  "=der:g5gyDx:asy=\0", "=der:g5gzDx:asy=\0",
  "=der:g5gzDy:asy=\0",
  //
  "=der:g5si00Dx:sym=\0", "=der:g5si00Dy:sym=\0", "=der:g5si00Dz:sym=\0",
  "=der:g5si0xDx:sym=\0", "=der:g5si0xDy:sym=\0", "=der:g5si0xDz:sym=\0",
  "=der:g5si0yDy:sym=\0", "=der:g5si0yDz:sym=\0",
  "=der:g5si0zDz:sym=\0",
  //
  "=der:g5six0D0:sym=\0", "=der:g5six0Dx:sym=\0", "=der:g5six0Dy:sym=\0", "=der:g5six0Dz:sym=\0",
  "=der:g5sixxDy:sym=\0", "=der:g5sixxDz:sym=\0",
  "=der:g5sixyDy:sym=\0", "=der:g5sixyDz:sym=\0",
  "=der:g5sixzDz:sym=\0",
  //
  "=der:g5siy0D0:sym=\0", "=der:g5siy0Dx:sym=\0", "=der:g5siy0Dy:sym=\0", "=der:g5siy0Dz:sym=\0",
  "=der:g5siyxDx:sym=\0", "=der:g5siyxDy:sym=\0", "=der:g5siyxDz:sym=\0",
  "=der:g5siyyDz:sym=\0",
  "=der:g5siyzDz:sym=\0",  
  //
  "=der:g5siz0D0:sym=\0", "=der:g5siz0Dx:sym=\0", "=der:g5siz0Dy:sym=\0", "=der:g5siz0Dz:sym=\0",
  "=der:g5sizxDx:sym=\0", "=der:g5sizxDy:sym=\0", "=der:g5sizxDz:sym=\0",
  "=der:g5sizyDy:sym=\0", "=der:g5sizyDz:sym=\0",
};

static enum {
  one, g5, g0, gx, gy, gz,
  g5g0, g5gx, g5gy, g5gz,
  g5si0x, g5si0y, g5si0z, g5sixy, g5sixz, g5siyz,
  noe_g0, noe_gx, noe_gy, noe_gz,
  der_g0_0, der_gx_0, der_gy_0, der_gz_0,
  der_gx_x, der_gy_x, der_gz_x,   
  der_gy_y, der_gz_y,
  der_gz_z,     
  der_g5g0_0_sym, der_g5gx_0_sym, der_g5gy_0_sym, der_g5gz_0_sym,
  der_g5gx_x_sym, der_g5gy_x_sym, der_g5gz_x_sym,   
  der_g5gy_y_sym, der_g5gz_y_sym,
  der_g5gz_z_sym,     
  der_g5gx_0_asy, der_g5gy_0_asy, der_g5gz_0_asy,
  der_g5gy_x_asy, der_g5gz_x_asy,   
  der_g5gz_y_asy,

  der_g5si00_x, der_g5si00_y, der_g5si00_z,
  der_g5si0x_x, der_g5si0x_y, der_g5si0x_z,
  der_g5si0y_y, der_g5si0y_z,
  der_g5si0z_z,

  der_g5six0_0, der_g5six0_x, der_g5six0_y, der_g5six0_z,
  der_g5sixx_y, der_g5sixx_z,
  der_g5sixy_y, der_g5sixy_z,
  der_g5sixz_z,

  der_g5siy0_0, der_g5siy0_x, der_g5siy0_y, der_g5siy0_z,
  der_g5siyx_x, der_g5siyx_y, der_g5siyx_z,
  der_g5siyy_z,
  der_g5siyz_z,  

  der_g5siz0_0, der_g5siz0_x, der_g5siz0_y, der_g5siz0_z,
  der_g5sizx_x, der_g5sizx_y, der_g5sizx_z,
  der_g5sizy_y, der_g5sizy_z,
} channels;

#endif /* _QHG_NN_THRP_DEFS_H */
