#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/sysinfo.h>
#include <stddef.h>
#include <unistd.h>
#include "thomas.h"

/*
* LMECA2660 Convection-Diffusion Equation
* by Thanh-Son TRAN 8116-12-00
*
*/

void initU(double* U_value,double sigma,double h,int N){
  double Q     = 1.0;
  for(int i = 0; i < N; i++){
    U_value[i] = Q/sqrt(M_PI*sigma*sigma) * exp(-((i-N/2)*(i-N/2)*h*h)/(sigma*sigma));
  }
  return;
}

void initExactU(double* U_value,double sigma,double h,double c,double L,double t,double nu,int N){
  // Note: Mathieu m'avait aider pour la logique de cette partie, qui m'a coûté par mal de temps
  //       pour rien. Du coups, ca logique doit sûrement transparaitre dans ce bout de code.
  double Q = 1.0;
  for(int i = 0; i < N; i++){
    double xtilde = ( i-N/2.0)*h - c*t + 1.0/2.0 * L;
    double xhat   = xtilde - floor( xtilde/L )*L;
    U_value[i] = Q/sqrt(M_PI*(sigma*sigma + 4.0*nu*t)) * exp( -( pow( xhat - L/2 ,2))/(sigma*sigma+4.0*nu*t));
  }
  return;
}

void cpyVector(double* A,double* B,int size){
    for(int i = 0; i< size; i++) B[i] = A[i];
    return;
}
double sumArray(double *A,int size){
  double inc = 0;
  for(int i = 0; i < size; i++){
    inc = inc + A[i];
  }
  return inc;
}
double sumArraySquare(double *A,int size){
  double inc = 0;
  for(int i = 0; i < size; i++){
    inc = inc + A[i]*A[i];
  }
  return inc;
}
void sumVector(double* A,double* B,double* C,int size){
    for(int i = 0; i< size; i++) C[i] = A[i] + B[i];
    return;
}
void diffVector(double* A,double* B,int size,double* C){
    for(int i = 0; i< size ;  i++) C[i] = A[i]-B[i];
    return;
}
void sum_Csum_Vector(double* A,double* B,double* C,int size,double inc){
    for(int i = 0; i< size; i++) A[i] = B[i]+C[i]*inc;
    return;
}
/*
* Solver of Finite Difference Convection (and) Diffusion Equation
*
*
* Note: Bon, il y a beaucoup de redodance, mais c'est pour le bien de la
*       compréhension de la matière, biensur il y a toujours moyen de
*       tout généraliser et faire cela bien proprement.
*/

void solverFDCEE2(double* U, double* du,double c,double h,double dt, int N){
  for(int m = 0; m < N; m++) du[m] =  dt*(-c)*(U[(m+1)%N] - U[(m+N-1)%N]) / (2.0*h);
  return;
}
void solverFDCEE4(double* U, double* du,double c,double h,double dt, int N){
  for(int m = 0; m < N; m++) du[m] =  dt*(-c)*(-1.0*U[(m+2)%N]+8.0*U[(m+1)%N]-8.0*U[(m-1+N)%N]+U[(m-2+N)%N])/(12.0*h);
  return;
}
void solverFDCES3(double* U, double* du,double c,double h,double dt, int N){
  for(int m = 0; m < N; m++) du[m] = dt*(-c)*(2.0*U[(m+1)%N]+3.0*U[(m)%N]-6.0*U[(m-1+N)%N]+U[(m-2+N)%N])/(6.0*h);
  return;
}
void solverFDCEI4(double* U, double* du,double c,double h,double dt, int N){
  double* q = calloc(N,sizeof(double));
  for(int m = 0; m < N ; m++){
    q[m] = 3.0/2.0*dt*(-c)*(U[(m+1)%N]-U[(m-1+N)%N])/(2.0*h) ;
  }
  solve_period_3diag(N,1,1.0/4.0,1.0/4.0,du,q);
  free(q);
  return;
}
void solverFDCEI6(double* U, double* du,double c,double h,double dt, int N){
  double* q = calloc(N,sizeof(double));
  for(int m = 0; m < N ; m++){
    q[m] = dt*(-c)*( 14.0/9.0*(U[(m+1)%N]-U[(m-1+N)%N])/(2.0*h) + 1.0/9.0*(U[(m+2)%N]-U[(m-2+N)%N])/(4.0*h) )  ;
  }
  solve_period_3diag(N,1,2.0/6.0,2.0/6.0,du,q);
  free(q);
  return;
}

void solverFDDEE2(double* U, double* d2u,double nu,double h,double dt, int N){
  for(int m = 0; m < N; m++) d2u[m] = dt*(nu)*(U[(m+1)%N]-2.0*U[(m)%N]+U[(m+N-1)%N])/(h*h);
  return;
}
