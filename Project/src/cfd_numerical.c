#include "cfd.h"


double scalar_psi_star_compute(struct _problem* Problem,int i,int j){
  return 0.5*( pow((*Problem).h,4.0) )/(2*pow((*Problem).h,2))*(
                 (*Problem).omega[i][j]
              + ((*Problem).psi[i+1][j] + (*Problem).psi[i-1][j]) /pow((*Problem).h,2)
              + ((*Problem).psi[i][j+1] + (*Problem).psi[i][j-1]) /pow((*Problem).h,2) );
}

double scalar_psi_compute(struct _problem* Problem,int i,int j){
  return (1-(*Problem).phi)*(*Problem).psi[i][j]+(*Problem).phi*scalar_psi_star_compute(*Problem,i,j);
}

void u_v_compute(struct _problem* Problem){
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
  for(int i; i < (*Problem).Nx ; i++){
    for(int j; j < (*Problem).Ny ; j++){

    }
  }
}
