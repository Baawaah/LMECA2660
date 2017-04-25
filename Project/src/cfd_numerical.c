#include "cfd.h"

void scalar_rhs(struct _problem* Problem, int i, int j, double du){
    double omega_loc = (*Problem).omega;
    double h_loc = (*Problem).h;
    double domdy,dom2dy2,domdx,dom2dx2 ;
    
    domdy = (omega_loc[i][j+1]-omega_loc[i][j-1])/(2*h_loc);
    dom2dy2 = (omega_loc[i][j+1]-2*omega_loc[i][j]+omega_loc[i][j-1])/h_loc^2;
    
    domdx = (omega_loc[i+1][j]-omega_loc[i-1][j])/(2*h_loc);
    dom2dx2 = (omega_loc[i+1][j]-2*omega_loc[i][j]+omega_loc[i-1][j])/h_loc^2;
    
    du = domdy + dom2dy2 + domdx + dom2dx2;
    
}

void first_iteration_omega(struct _problem* Problem, double du_old){
    for(int i=0; i<(*Problem).Nx; i++){
        for(int j = 0; j < (*Problem).Ny ; j++){
            du_old=scalar_rhs((*Problem),i,j,du);
            (*Problem).omega[i][j]=(*Problem).omega[i][j]+du_old*(*Problem).dt;
        }
    }
}

voir integration_omega(struct _problem* Problem){
    double du,du_old;
    first_iteration_omega(*Problem,du_old);
    
    for(int i=0; i<(*Problem).Nx; i++){
        for(int j = 0; j < (*Problem).Ny ; j++){
            du=scalar_rhs((*Problem),i,j,du);
            (*Problem).omega[i][j]=(*Problem).omega[i][j]+dt/2*(3*du-du_old);
        }
    }
    
}
double scalar_psi_star_compute(struct _problem* Problem,int i,int j){
  return 0.5*( pow((*Problem).h,4.0) )/(2*pow((*Problem).h,2))*(
                 (*Problem).omega[i][j]
              + ((*Problem).psi[i+1][j] + (*Problem).psi[i-1][j]) /pow((*Problem).h,2)
              + ((*Problem).psi[i][j+1] + (*Problem).psi[i][j-1]) /pow((*Problem).h,2) );
}
double scalar_psi_compute(struct _problem* Problem,int i,int j){
  return (1-(*Problem).phi)*(*Problem).psi[i][j]+(*Problem).phi*scalar_psi_star_compute(Problem,i,j);
}
double scalar__psi_r_compute(struct _problem* Problem,int i, int j){
  return ((*Problem).psi[i+1][j] - 4*(*Problem).psi[i][j] + (*Problem).psi[i-1][j] +(*Problem).psi[i][j+1] + (*Problem).psi[i][j-1])/(pow((*Problem).h,2)) + (*Problem).omega[i][j];
}

void inner_u_v_compute(struct _problem* Problem){
  for(int i; i < (*Problem).Nx ; i++){
     for(int j; j < (*Problem).Ny ; j++){
       // by centred finite differences
       (*Problem).u    [i][j] = ((*Problem).psi[i+1][j]-(*Problem).psi[i-1][j])/2*(*Problem).h;
       (*Problem).v    [i][j] = ((*Problem).psi[i][j+1]-(*Problem).psi[i][j-1])/2*(*Problem).h;
     }
  }
}

void boundary_psi_update(struct _problem* Problem, double (*Q)(double) ){
// Must be change in function of the actual domain
// Need a Q(t) function

// Upper boundary
  for(int i = 0; i < (*Problem).Nx ; i++ ) (*Problem).psi[i][0] = Q((*Problem).t);

// Down boundary - Left - Right
  for(int i = 0; i < (*Problem).Nx ; i++ ){
    if   (i < (*Problem).NLs ){ (*Problem).psi[i][(*Problem).NHs-1] = 0; }
    else                        (*Problem).psi[i][(*Problem).Ny -1] = 0;
  }
// Down boundary Side
  for(int i = (*Problem).NHs; i < (*Problem).Ny; i++ ) (*Problem).psi[(*Problem).NLs-1][i] = 0;

// Inflow Boundary - Natural Condition
  for(int i = 1; i < (*Problem).NHs-1 ; i++) (*Problem).psi[0][i]               = (*Problem).psi[1][i];

// Outflow Boundary - Natural Condition
  for(int i = 1; i < (*Problem).Ny -1 ; i++) (*Problem).psi[(*Problem).Nx-1][i] = (*Problem).psi[(*Problem).Nx-2][i];

}

void boundary_omega_update(struct _problem* Problem){
// Must be change in function of the actual domain

// Upper boundary
  for(int i = 0; i < (*Problem).Nx ; i++ ) (*Problem).omega[i][0] = -3.0/((*Problem).h*(*Problem).h) * (*Problem).psi[i][1] - 0.5*(*Problem).omega[i][1];

// Down boundary - Left - Right
  for(int i = 0; i < (*Problem).Nx ; i++ ){
    if   (i < (*Problem).NLs ){ (*Problem).omega[i][(*Problem).NHs-1] = -3.0/((*Problem).h*(*Problem).h) * (*Problem).psi[i][(*Problem).NHs-2] - 0.5*(*Problem).omega[i][(*Problem).NHs-1]; }
    else                        (*Problem).omega[i][(*Problem).Ny -1] = -3.0/((*Problem).h*(*Problem).h) * (*Problem).psi[i][(*Problem).Ny -2] - 0.5*(*Problem).omega[i][(*Problem).Ny -1];
  }
// Down boundary Side
  for(int i = (*Problem).NHs; i < (*Problem).Ny; i++ ) (*Problem).omega[(*Problem).NLs-1][i] = -3.0/((*Problem).h*(*Problem).h) * (*Problem).psi[(*Problem).NLs][i] - 0.5*(*Problem).omega[(*Problem).NLs][i];

// Corner boundary
  (*Problem).omega[(*Problem).NLs-1][(*Problem).NHs-1] = - 3.0 / (2*(*Problem).h*(*Problem).h) * (*Problem).psi[(*Problem).NLs-2][(*Problem).NHs-2] - 0.5* (*Problem).omega[(*Problem).NLs-2][(*Problem).NHs-2];
                                       // pas plutot ? = - 3.0 / (2*(*Problem).h*(*Problem).h) * (*Problem).psi[(*Problem).NLs][(*Problem).Ny   ] - 0.5* (*Problem).omega[(*Problem).NLs][(*Problem).NHs];
//(*Problem).omega[(*Problem).NLs-1][(*Problem).Ny -1] = - 3.0 / (2*(*Problem).h*(*Problem).h) * (*Problem).psi[(*Problem).NLs][(*Problem).Ny -2] - 0.5* (*Problem).omega[(*Problem).NLs][(*Problem).Ny -2];
                                       // pas plutot ? = - 3.0 / (2*(*Problem).h*(*Problem).h) * (*Problem).psi[(*Problem).NLs][(*Problem).Ny   ] - 0.5* (*Problem).omega[(*Problem).NLs][(*Problem).Ny   ];

// Inflow Boundary - Natural Condition
  for(int i = 1; i < (*Problem).NHs-1 ; i++) (*Problem).omega[0][i]               = (*Problem).omega[1][i];

// Outflow Boundary - Natural Condition
  for(int i = 1; i < (*Problem).Ny -1 ; i++) (*Problem).omega[(*Problem).Nx-1][i] = (*Problem).omega[(*Problem).Nx-2][i];

}

void inner_psi_star_update(struct _problem* Problem){
  for(int i = 1; i < (*Problem).Nx -1; i++){
    for(int j = 1; j < (*Problem).Ny -1; j++){
      (*Problem).psi[i][j] = scalar_psi_compute(Problem,i,j);
    }
  }
}

double inner_psi_error_compute(struct _problem* Problem){
  double e_error = 0;
  for(int i=1; i < (*Problem).Nx -1; i++){
    for(int j=1; j < (*Problem).Ny -1; j++){
      double R = scalar__psi_r_compute(Problem,i,j);
      e_error = (*Problem).H/(*Problem).Q0 * sqrt(R*R*(*Problem).h*(*Problem).h) + e_error;
    }
  }
  return e_error;
}

void inner_psi_interator(struct _problem* Problem){
  while( (*Problem).e_max < inner_psi_error_compute(Problem) ){
    inner_psi_star_update(Problem);
  }
}
