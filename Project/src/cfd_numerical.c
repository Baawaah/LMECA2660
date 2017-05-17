#include "cfd.h"
#include "stdio.h"
#include <stddef.h>

double scalar_psi_star_compute(struct _problem* Problem,int i,int j){
  return 0.25*(
                 (*Problem).omega[i][j] * ((*Problem).h*(*Problem).h)
              + ((*Problem).psi[i+1][j] + (*Problem).psi[i-1][j])
              + ((*Problem).psi[i][j+1] + (*Problem).psi[i][j-1])    );
}

double scalar_psi_compute(struct _problem* Problem,int i,int j){
  return    (1.0-(*Problem).phi) * (*Problem).psi[i][j]
          +      (*Problem).phi  * scalar_psi_star_compute(Problem,i,j);
}

double scalar_psi_r_compute(struct _problem* Problem,int i, int j){
  return   (*Problem).omega[i][j] +
          (  (*Problem).psi[i+1][j]   + (*Problem).psi[i-1][j]
            - 4*(*Problem).psi[i][j]
            +(*Problem).psi[i]  [j+1] + (*Problem).psi[i]  [j-1]
          )/pow((*Problem).h,2);
}

double scalar_u_v_poiseuille     (struct _problem* Problem,double y){
  return 6.0*functionQ(Problem)/((*Problem).Hc*(*Problem).Hc) *  ( y            - (y*y)/(*Problem).Hc              );
}
double scalar_u_v_poiseuille_dy  (struct _problem* Problem,double y){
  return -6.0*functionQ(Problem)/((*Problem).Hc*(*Problem).Hc)  * ( 1.0          - (2.0*y)/(*Problem).Hc            );
}
double scalar_u_v_poiseuille_int (struct _problem* Problem,double y){
  return 6.0*functionQ(Problem)/((*Problem).Hc*(*Problem).Hc)  * ((1.0/2.0)*y*y - (1.0/3.0)*(y*y*y)/(*Problem).Hc) + functionQ(Problem);
}


double scalar_u_v_poiseuille_eta(struct _problem* Problem,double eta){
  return 1.5*functionQ(Problem)*(1.0 - pow(eta,2));
}

double scalar_u_v_poiseuille_eta_dy(struct _problem* Problem,double eta){
  return -3.0*functionQ(Problem)*(eta);
}

double scalar_u_v_poiseuille_eta_int(struct _problem* Problem,double eta){
  return 3.0/2.0*functionQ(Problem)*(eta - eta*eta*eta * 1.0/3.0) + functionQ(Problem);
}

void inner_u_v_compute(struct _problem* Problem){
  for(int i= 1 ; i < (*Problem).Nx-1; i++){
     for(int j=(*Problem).imap[i]+1; j < (*Problem).Ny-1; j++){
       // by centred finite differences
       (*Problem).u    [i][j] =     ((*Problem).psi[i][j+1]-(*Problem).psi[i][j-1])/(2.0*(*Problem).h);
       (*Problem).v    [i][j] =    -((*Problem).psi[i+1][j]-(*Problem).psi[i-1][j])/(2.0*(*Problem).h);
     }
  }
}

void u_v_stag(struct _problem* Problem){
  for(int i= 1 ; i < (*Problem).Nx_p-1; i++){
     for(int j=(*Problem).imap_p[i]+1; j < (*Problem).Ny_p-1; j++){
        (*Problem).u    [i][j] =     ((*Problem).psi[i][j+1]-(*Problem).psi[i][j])/(*Problem).h;
        (*Problem).v    [i][j] =    -((*Problem).psi[i+1][j]-(*Problem).psi[i][j])/(*Problem).h;
     }
  }
}

void boundary_psi_update(struct _problem* Problem, double (*Q)(struct _problem*)){
// Upper boundary
  for(int i = 0; i < (*Problem).Nx ; i++ ) (*Problem).psi[i][(*Problem).Ny-1] = Q(Problem);

// Down boundary - Left - Right
  for(int i = 0; i < (*Problem).Nx ; i++ ){
    (*Problem).psi[i][(*Problem).imap[i]] = 0.0;
  }
// Down boundary Side
  if( (*Problem).Ls != (*Problem).L && (*Problem).Ls != 0.0 )
    for(int j = 0; j < (*Problem).NHs; j++ ) (*Problem).psi[(*Problem).NLs-1][j] = 0.0;

}

void boundary_omega_update(struct _problem* Problem){
// Upper boundary
  for(int i = 0; i < (*Problem).Nx ; i++ ) (*Problem).omega[i][(*Problem).Ny-1   ] = -3.0/((*Problem).h*(*Problem).h) * ((*Problem).psi[i][(*Problem).Ny-2    ] - functionQ(Problem)) - 0.5*(*Problem).omega[i][(*Problem).Ny-2     ];

// Down boundary - Left - Right
  for(int i = 0; i < (*Problem).Nx ; i++ ) (*Problem).omega[i][(*Problem).imap[i]] = -3.0/((*Problem).h*(*Problem).h) * (*Problem).psi[i][(*Problem).imap[i]+1] - 0.5*(*Problem).omega[i][(*Problem).imap[i]+1];

// Down boundary Side
  if( (*Problem).Ls != (*Problem).L && (*Problem).Ls != 0.0 )
  for(int j = 0; j < (*Problem).NHs; j++ ) (*Problem).omega[(*Problem).NLs-1][j] = -3.0/((*Problem).h*(*Problem).h) * (*Problem).psi[(*Problem).NLs+1][j] - 0.5*(*Problem).omega[(*Problem).NLs+1][j];
// Corner boundary
  if( (*Problem).Ls != (*Problem).L && (*Problem).Ls != 0.0 )
  (*Problem).omega[(*Problem).NLs-1][(*Problem).NHs-1] = - 3.0 / (2*(*Problem).h*(*Problem).h) * (*Problem).psi[(*Problem).NLs][(*Problem).NHs] - 0.5* (*Problem).omega[(*Problem).NLs][(*Problem).NHs];

}
void boundary_omega_dwdx_update(struct _problem* Problem){
  // Inflow Boundary - Natural Condition
  for(int j = (*Problem).imap[0]; j < (*Problem).Ny-1 ; j++) (*Problem).omega[0][j] = (*Problem).omega[1][j];
  // Outflow Boundary - Natural Condition
  for(int i = (*Problem).imap[(*Problem).Nx-1]; i < (*Problem).Ny-1 ; i++) (*Problem).omega[(*Problem).Nx-1][i] = (*Problem).omega[(*Problem).Nx-2][i];
}

void inner_psi_update(struct _problem* Problem){

  for(int i = 1; i < (*Problem).Nx-1; i++){
    for(int j = (*Problem).imap[i]+1; j < (*Problem).Ny-1; j++){
      (*Problem).psi[i][j] = scalar_psi_compute(Problem,i,j);
    }
  }/*
  for(int i = (*Problem).Nx-2; i > 0 ; i--){
   for(int j = 1; j < (*Problem).imax_map[i]-1; j++){
      (*Problem).psi[i][j] = scalar_psi_compute(Problem,i,j);
    }
  }*/
}

double inner_psi_error_compute(struct _problem* Problem){
  double e_error = 0.0;
  double square = 0.0;
  for(int i=1; i < (*Problem).Nx -1; i++){
    for(int j= (*Problem).imap[i]+1; j < (*Problem).Ny-1; j++){
      (*Problem).R_res[i][j] = scalar_psi_r_compute(Problem,i,j);
      square = (*Problem).R_res[i][j]*(*Problem).R_res[i][j] + square;
    }
  }
  e_error = fabs((*Problem).H*(*Problem).H/(*Problem).Q0*(*Problem).h*sqrt(1.0/((*Problem).L*(*Problem).H) *square));
  return e_error;
}

void poisson_inner_psi_iterator(struct _problem* Problem){
  int n_iter = 0;
  int iter = 0 ;
  int iter_max = 5000;
  double error = (*Problem).tol+1;

  while( error>(*Problem).tol && iter < iter_max){
    n_iter++;
    inner_psi_update(Problem);
    // Inflow Boundary - Natural Condition
    for(int j = (*Problem).imap[0]; j < (*Problem).Ny-1 ; j++) (*Problem).psi[0][j]               = (*Problem).psi[1][j];
    // Outflow Boundary - Natural Condition
    for(int i = (*Problem).imap[(*Problem).Nx-1]; i < (*Problem).Ny -1 ; i++) (*Problem).psi[(*Problem).Nx-1][i] = (*Problem).psi[(*Problem).Nx-2][i];

    error = inner_psi_error_compute(Problem);
    iter++;
  }
  if(iter >= iter_max){fprintf(stderr, "[DEADSTOP] Maximum Iteration Reached Current Error: %f\n",error); deadstop_exit(Problem);}
}

int diagnose_check(struct _problem* Problem,int i,int j,int ktime){
  int check = 0;
  double Re_h       = (*Problem).h*(fabs((*Problem).u[i][j]) + fabs((*Problem).v[i][j]))/(*Problem).nu;
  if(Re_h > (*Problem).Re_h[ktime]) (*Problem).Re_h[ktime] = Re_h;
  double Re_h_omega = (*Problem).h*(*Problem).h*fabs((*Problem).omega[i][j])/(*Problem).nu;
  if(Re_h_omega > (*Problem).Re_h_omega[ktime]) (*Problem).Re_h_omega[ktime] = Re_h_omega;
  double Beta       = (*Problem).h*(fabs((*Problem).u[i][j]) + fabs((*Problem).v[i][j]))*(*Problem).dt /(*Problem).h;
  if(Beta > (*Problem).Beta_CFL[ktime]) (*Problem).Beta_CFL[ktime] = Beta;
  if( Re_h       >=  5 )                  check += 1;
  if( Re_h_omega >= 20 )                  check += 2;
  if( Beta       >=  (*Problem).CFL    )  check += 4;
  return check;
}
