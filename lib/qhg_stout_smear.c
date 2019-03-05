#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
#include <qhg_defs.h>
#include <qhg_su3_ops.h>
#include <qhg_su3_project.h>
#include <qhg_gauge_field.h>
#include <qhg_xchange_gauge.h>
#include <math.h>

#define G(v,mu) (((v)*ND + mu)*NC*NC)

void qhg_su3_print(_Complex double *in)
{
  printf("\n");    

  for(int c0=0; c0<NC; c0++){
    for(int c1=0; c1<NC; c1++)
      printf("%f%+fi,",creal(in[c0*NC+c1]),cimag(in[c0*NC+c1]));

    printf("\n");    
  }
  printf("\n");    

  return;
}

/* needed function for the exponentiation */
double qhg_xi0(double w)
{
  double tmp,tmp1;
  if(abs(w) < 0.05)
    {
      tmp=1.0-w*w/42.0;
      tmp1=1.0-w*w/20.0*tmp;
      tmp=1.0-w*w/6.0*tmp1;
      return tmp;
    }
  else
    return sin(w)/w;
}

/* analytic exponentiation: Peardon & Morningstar (arXiv:hep-lat/0311018) */
void qhg_exp(_Complex double *M)
{
  _Complex double c0,c1;
  double c0max,u,w,theta;
  _Complex double f0,f1,f2;
  _Complex double h0,h1,h2;

  _Complex double M2[NC*NC],M3[NC*NC],one[NC*NC];

  _Complex double C,C1;
  _Complex double tmp,tmp1;
  double factor;
  int sign=0;

  su3_mul_uu(M2,M,M); //M2=M*M
  su3_mul_uu(M3,M2,M); //M3=M*M*M

  c0=su3_linalg_trace_u(M3);
  c1=su3_linalg_trace_u(M2);

  c0 = c0/3.0;
  c1 = c1/2.0;
  
  if(creal(c0)<0){
    sign=1;
    c0 *= -1.0;
  }
  c0max=2.0*pow(c1/3.0,1.5);

  theta=acos(c0/c0max);

  u=sqrt(c1/3.0)*cos(theta/3.0);
  w=sqrt(c1)*sin(theta/3.0);

  C=cos(2*u) + I*sin(2*u); //exp(2iu)
  C1=cos(u) - I*sin(u);   //exp(-iu)

  tmp=u*u-w*w;
  tmp1=8*u*u*cos(w) + I*2*u*(3*u*u+w*w)*qhg_xi0(w);

  h0=tmp*C+C1*tmp1; //final h0,eq.(30)

  tmp=2*u;
  tmp1=2*u*cos(w) - I*(3*u*u-w*w)*qhg_xi0(w);
  h1=tmp*C-C1*tmp1; //final h1, eq.(31)

  tmp=cos(w) + I*3*u*qhg_xi0(w);
  h2=C-C1*tmp; //final h2, eq.(32)

  factor=9*u*u-w*w;  //eq.(29)

  if(sign==0) //i.e. creal(c0)>0
    {
      f0=h0/factor; // Eq.(29)
      f1=h1/factor;
      f2=h2/factor;
    }
  else
    {
      f0=conj(h0/factor); // Eq.(34)
      f1=conj(-h1/factor);
      f2=conj(h2/factor);
    }

  _Complex double f0_arr[NC] = {f0,f0,f0};
  su3_linalg_diag(one, f0_arr); //one=f0*I
  su3_linalg_au(f1,M);          //f1*M and store in M
  su3_linalg_au(f2,M2);         //f2*M^2

  su3_linalg_upeqv(M,one); //M=M+one
  su3_linalg_upeqv(M,M2); //M=M+M2=f0*I+f1*Q+f2*Q^2 -> exp(iQ), eq.(19)

  return;
}

void
qhg_stout_smear_iter(qhg_gauge_field out, qhg_gauge_field in, double omega)
{
#ifdef QHG_OMP
#pragma omp parallel
  {
#endif
    
  unsigned long int vol = in.lat->vol;
  unsigned long int lvol = in.lat->lvol;
  unsigned long int **nn = in.lat->nn;

#ifdef QHG_OMP
#pragma omp single
  {
#endif  

  qhg_xchange_gauge(in);

#ifdef QHG_OMP
  }
#endif  

  _Complex double *U = in.field;
  _Complex double *V = out.field;  
  _Complex double *u0, *u1, *u2;
  _Complex double staple[NC*NC], w[NC*NC], u[NC*NC];

  for(int mu=0; mu<ND; mu++)
#ifdef QHG_OMP
#pragma omp for
#endif
    for(unsigned long int v00=0; v00<lvol; v00++) {
      su3_linalg_zero(staple);

      /* create the staple term for this mu direction */
      for(int nu=(mu==0?1:0); nu<ND; nu==mu-1?nu+=2:nu++){
      	unsigned long int vp0 = nn[mu][v00];
      	unsigned long int v0p = nn[nu][v00];
      	unsigned long int v0m = nn[nu+ND][v00];	
      	unsigned long int vpm = nn[nu+ND][nn[mu][v00]];
	
      	/* Fwd staple */
      	u0 = &U[G(v00, nu)];
      	u1 = &U[G(v0p, mu)];
      	u2 = &U[G(vp0, nu)];
      	su3_mul_uu(u, u0, u1);
      	su3_mul_ud(w, u, u2);
	
      	su3_linalg_upeqv(staple, w);

      	/* Bwd staple */
      	u0 = &U[G(v0m, nu)];
      	u1 = &U[G(v0m, mu)];
      	u2 = &U[G(vpm, nu)];
      	su3_mul_du(u, u0, u1);
      	su3_mul_uu(w, u, u2);
	
      	su3_linalg_upeqv(staple, w);
      }
      u0 = &U[G(v00, mu)];
      u1 = &V[G(v00, mu)];

      /* Omega = staple*U^+ */

      su3_mul_ud(w, staple, u0);

      /* Omega^+ - Omega */
      su3_linalg_ueqvd(u, w);
      su3_linalg_umeqv(u, w);
      
      /* calculate trace term of Q 
	 1/NC * tr(Omega^+ - Omega) */

      _Complex double trace = su3_linalg_trace_u(u);
      _Complex double factor[NC] = {trace/NC, trace/NC, trace/NC};
      su3_linalg_diag(u2, factor);

      /* calculate Q 
	 Q = i*omega/2*[(Omega^+ - Omega) - 1/NC*tr(Omega^+ - Omega)] */

      su3_linalg_umeqv(u, u2);
      su3_linalg_au(I*0.5*omega, u);

      /* calculate exp(iQ) */

      qhg_exp( u );

      /* Multiply exp(iQ) by U */

      su3_mul_uu(u1, u, u0);

      if(in.lat->comms->proc_id==0 && mu==0 && v00==0){
	printf("U^(n+1):\n");
	qhg_su3_print(u1);
      }

    }
#ifdef QHG_OMP
  }
#endif
  return;
}

void
qhg_stout_smear(qhg_gauge_field out, qhg_gauge_field in, double omega, int niter)
{
  qhg_gauge_field aux[2];
  aux[0] = qhg_gauge_field_init(in.lat);
  aux[1] = qhg_gauge_field_init(in.lat);

  qhg_gauge_field_copy(aux[0], in);
  qhg_gauge_field_copy(aux[1], in);
  for(int i=0; i<niter; i++)
    qhg_stout_smear_iter(aux[(i+1) % 2], aux[i % 2], omega);
  qhg_gauge_field_copy(out, aux[niter % 2]);
  
  qhg_gauge_field_finalize(aux[0]);
  qhg_gauge_field_finalize(aux[1]);  
  return;
}
