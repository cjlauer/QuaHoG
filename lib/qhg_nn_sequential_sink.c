#include <qhg_types.h>
#include <qhg_idx.h>
#include <qhg_prop_ops.h>
#include <qhg_prop_proj.h>
#include <qhg_prop_gammas.h>
#include <math.h>

void
qhg_nn_sequential_sink_u(qhg_spinor_field out[], qhg_spinor_field in_up[], qhg_spinor_field in_dn[],
			 int t0, qhg_thrp_nn_sink_params thrp_params)
{
  qhg_lattice *lat = in_up[0].lat;
  int *proc_coords = lat->comms->proc_coords;
  int *dims = lat->dims;
  int lt = lat->ldims[0];  
  unsigned long int lvol = lat->lvol;
  unsigned long int lv3 = lat->lv3;

  /* Only loop along appropriate sink time-slice */
  int dt = thrp_params.dt;
  enum projector proj = thrp_params.proj;
  int ts = (t0+dt) % dims[0];
  int t_rank = ts / lt;
  int ts_loc = ts % lt;
  unsigned long int v0, v1;
  if(proc_coords[0] == t_rank) {
    v0 = lv3*(unsigned long int)ts_loc;
    v1 = lv3*(unsigned long int)(ts_loc+1);
  } else {
    v0 = 0; /* Do nothing */
    v1 = 0;
  }

  _Complex double zero[NC*NS][NC*NS];
  _Complex double one[NC*NS][NC*NS];
  prop_zero(zero);
  prop_eye(one);
  for(unsigned long int v=0; v<lvol; v++)
    prop_store(out, v, zero);
  
  _Complex double P[NS*NC][NS*NC];
  switch(proj) {
  case P0:
    prop_proj_TP0T_G(P, one);
    break;
  case P3:
    prop_proj_TP3T_G(P, one);
    break;
  case P4:
    prop_proj_TP4T_G(P, one);
    break;
  case P5:
    prop_proj_TP5T_G(P, one);
    break;
  case P6:
    prop_proj_TP6T_G(P, one);
    break;
  }
  
  for(unsigned long int v=v0; v<v1; v++) {
    /* U is in_up[], D is in_dn[], A0, B0, A1, B1, X are auxiliary, Y is
       out[]. */
    _Complex double U[NS*NC][NS*NC];
    _Complex double D[NS*NC][NS*NC];    
    _Complex double X[NS*NC][NS*NC];    
    _Complex double Y[NS*NC][NS*NC];    
    _Complex double A0[NS*NC][NS*NC];
    _Complex double B0[NS*NC][NS*NC];    
    _Complex double A1[NS*NC][NS*NC];
    _Complex double B1[NS*NC][NS*NC];    
    prop_load(U, in_up, v);
    prop_load(D, in_dn, v);
    prop_zero(A0);
    
    _Complex double Pu[NS*NC][NS*NC];
    _Complex double uP[NS*NC][NS*NC];
    _Complex double Cg5dCg5[NS*NC][NS*NC];
    switch(proj) {
    case P0:
      prop_proj_TP0T_G(Pu, U);
      prop_proj_G_TP0T(uP, U);
      break;
    case P3:
      prop_proj_TP3T_G(Pu, U);
      prop_proj_G_TP3T(uP, U);
      break;
    case P4:
      prop_proj_TP4T_G(Pu, U);
      prop_proj_G_TP4T(uP, U);
      break;
    case P5:
      prop_proj_TP5T_G(Pu, U);
      prop_proj_G_TP5T(uP, U);
      break;
    case P6:
      prop_proj_TP6T_G(Pu, U);
      prop_proj_G_TP6T(uP, U);
      break;
    }
    
    prop_Cg5_G(X, D);
    prop_G_Cg5(Cg5dCg5, X);
    
    prop_contract_02(X, Cg5dCg5, U);
    prop_contract_23(A1, Cg5dCg5, uP);
    prop_contract_02(B0, uP, Cg5dCg5);
    prop_contract_13(B1, Cg5dCg5, Pu);
    for(int mu=0; mu<NS; mu++)
      for(int nu=0; nu<NS; nu++)
	for(int c0=0; c0<NC; c0++)
	  for(int c1=0; c1<NC; c1++)
	    for(int ku=0; ku<NS; ku++) {
	      A0[CS(mu, c0)][CS(nu, c1)] += X[CS(ku, c0)][CS(ku, c1)]*P[CS(nu, 0)][CS(mu, 0)];
	    }
    
    for(int cs0=0; cs0<NS*NC; cs0++)
      for(int cs1=0; cs1<NS*NC; cs1++)
	X[cs0][cs1] = conj(A0[cs0][cs1] + A1[cs0][cs1] + B0[cs0][cs1] + B1[cs0][cs1]);
    
    prop_g5_G(Y, X);
    prop_store(out, v, Y);
  }
  return;
}

void
qhg_nn_sequential_sink_d(qhg_spinor_field out[], qhg_spinor_field in[], int t0, qhg_thrp_nn_sink_params thrp_params)
{
  qhg_lattice *lat = in[0].lat;
  int *proc_coords = lat->comms->proc_coords;
  int *dims = lat->dims;
  int lt = lat->ldims[0];  
  unsigned long int lvol = lat->lvol;
  unsigned long int lv3 = lat->lv3;

  /* Only loop along appropriate sink time-slice */
  int dt = thrp_params.dt;
  enum projector proj = thrp_params.proj;
  int ts = (t0+dt) % dims[0];
  int t_rank = ts / lt;
  int ts_loc = ts % lt;
  unsigned long int v0, v1;
  if(proc_coords[0] == t_rank) {
    v0 = lv3*(unsigned long int)ts_loc;
    v1 = lv3*(unsigned long int)(ts_loc+1);
  } else {
    v0 = 0; /* Do nothing */
    v1 = 0;
  }

  _Complex double zero[NC*NS][NC*NS];
  prop_zero(zero);
  for(unsigned long int v=0; v<lvol; v++)
    prop_store(out, v, zero);
  
  for(unsigned long int v=v0; v<v1; v++) {
    /* X is in[], Y is out[], A, B are auxiliary */
    _Complex double X[NS*NC][NS*NC];
    _Complex double Y[NS*NC][NS*NC];
    _Complex double A[NS*NC][NS*NC];
    _Complex double B[NS*NC][NS*NC];    
    prop_load(X, in, v);
    _Complex double PxCg5[NS*NC][NS*NC];
    _Complex double Cg5xCg5[NS*NC][NS*NC];
    _Complex double Cg5x[NS*NC][NS*NC];
    _Complex double Px[NS*NC][NS*NC];
    switch(proj) {
    case P0:
      prop_proj_TP0T_G(Px, X);
      break;
    case P3:
      prop_proj_TP3T_G(Px, X);
      break;
    case P4:
      prop_proj_TP4T_G(Px, X);
      break;
    case P5:
      prop_proj_TP5T_G(Px, X);
      break;
    case P6:
      prop_proj_TP6T_G(Px, X);
      break;
    }
    
    prop_Cg5_G(Cg5x, X);
    prop_G_Cg5(PxCg5, Px);
    prop_G_Cg5(Cg5xCg5, Cg5x);
    
    prop_contract_12(A, Cg5x, PxCg5);
    prop_contract_01(B, Px, Cg5xCg5);

    for(int cs0=0; cs0<NS*NC; cs0++)
      for(int cs1=0; cs1<NS*NC; cs1++)
	Y[cs0][cs1] = conj(A[cs0][cs1] + B[cs0][cs1]);
    
    prop_g5_G(X, Y);
    prop_store(out, v, X);
  }
  return;
}

void
qhg_mesons_sequential_sink(qhg_spinor_field out[], qhg_spinor_field in[], int t0, qhg_thrp_nn_sink_params thrp_params)
{
  qhg_lattice *lat = in[0].lat;
  int *proc_coords = lat->comms->proc_coords;
  int *dims = lat->dims;
  int lt = lat->ldims[0];  
  unsigned long int lvol = lat->lvol;
  unsigned long int lv3 = lat->lv3;

  /* Only loop along appropriate sink time-slice */
  int dt = thrp_params.dt;
  int ts = (t0+dt) % dims[0];
  int t_rank = ts / lt;
  int ts_loc = ts % lt;
  unsigned long int v0, v1;
  if(proc_coords[0] == t_rank) {
    v0 = lv3*(unsigned long int)ts_loc;
    v1 = lv3*(unsigned long int)(ts_loc+1);
  } else {
    v0 = 0; /* Do nothing */
    v1 = 0;
  }

  _Complex double zero[NC*NS][NC*NS];
  prop_zero(zero);
  for(unsigned long int v=0; v<lvol; v++)
    prop_store(out, v, zero);
  
  for(unsigned long int v=v0; v<v1; v++) {
    /* X is in[], Y is out[] */
    _Complex double X[NS*NC][NS*NC];
    _Complex double Y[NS*NC][NS*NC];
    prop_load(X, in, v);
    _Complex double g5x[NS*NC][NS*NC];
    
    prop_g5_G(Y, X);
    prop_store(out, v, Y);
  }
  return;
}

void
qhg_phase_sequential_sink(qhg_spinor_field out[], qhg_spinor_field in[], int mom[3], int source_coords[4], double sign)
{
  qhg_lattice *lat = out[0].lat;
  int *proc_coords = lat->comms->proc_coords;
  int *dims = lat->dims;
  int *ldims = lat->ldims;
  unsigned long int lvol = lat->lvol;

  double shift[3];
  for(int d=0; d<3; d++)
    shift[d] = dims[d+1] + ldims[d+1]*proc_coords[d+1] - source_coords[d+1];

  for(unsigned long int v=0; v<lvol; v++) {
    int co[4] = CO(v,ldims);
    _Complex double coeff = 0;
    double ph = 0;
    for(int d=0; d<3; d++)
      ph += (mom[d]*(shift[d]+co[d+1]))/dims[d+1];
    ph *= sign*2*M_PI;
    coeff = cos(ph) - _Complex_I * sin(ph);
    for(int cs0=0; cs0<NS*NC; cs0++)
      for(int cs1=0; cs1<NS*NC; cs1++) {
        out[cs0].field[VSC(v, cs1)] = in[cs0].field[VSC(v, cs1)]*coeff;
      }
  }
}
