#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/sysinfo.h>
#include <stddef.h>
#include <unistd.h>
#include "thomas.h"

void initU(double* U_value,double h,int N){
  double Q     = 1.0;
  double sigma = (2.0)*h;
  for(int i = 0; i < N; i++){
    U_value[i] = Q/sqrt(M_PI*sigma*sigma) * exp(-((i-N/2)*(i-N/2))/(sigma*sigma));
  }
  //U_value[N/2] = 10.0;
  return;
}

void cpyVector(double* A,double* B,int size){
    for(int i = 0; i< size; i++) B[i] = A[i];
    return;
}
void sumVector(double* A,double* B,int size,double inc){
    for(int i = 0; i< size; i++) B[i] = A[i]+inc;
    return;
}
void sum_Csum_Vector(double* A,double* B,double* C,int size,double inc){
    for(int i = 0; i< size; i++) A[i] = B[i]+C[i]*inc;
    return;
}

double FDEE2(double* U,double h,int N,int m){
//  Finite Difference O(h^2)
return (U[(m+1)%N] - U[(m+N-1)%N]) / (2.0*h);
}
double FDEE4(double* U,double h,int N,int m){
//  Finite Difference O(h^2)
return (U[(m+2)%N]+8*U[(m+1)%N]-8*U[(m-1+N)%N]+U[(m-1+N)%N])/(12*h);
}
void solverFDEE2(double* U, double* du,double c,double h,double dt, int N){
  for(int m = 0; m < N; m++) du[m] =  dt*(-c)*(U[(m+1)%N] - U[(m+N-1)%N]) / (2.0*h);
  return;
}
void solverFDEE4(double* U, double* du,double c,double h,double dt, int N){
  for(int m = 0; m < N; m++) du[m] =  dt*(-c)*(U[(m+2)%N]+8.0*U[(m+1)%N]-8.0*U[(m-1+N)%N]+U[(m-2+N)%N])/(12.0*h);
  return;
}
void solverFDES3(double* U, double* du,double c,double h,double dt, int N){
  for(int m = 0; m < N; m++) du[m] = dt*(-c)*(2.0*U[(m+1)%N]+3.0*U[(m)%N]-6.0*U[(m-1+N)%N]+U[(m-2+N)%N])/(12.0*h);
  return;
}

int main(int argc,char* argv[]){
// Init Value
  double c = 1.0;
  int N = 128;
  int Ntime = 1000;
  double h = 0.2;
  double dt = 0.2;
  double* U     = calloc(N    ,sizeof(double));
  //double* t     = calloc(N,sizeof(double));
  //double* dudx  = calloc(N    ,sizeof(double));
  double* du  = calloc(N ,sizeof(double));

  double* Us = calloc(N , sizeof(double));
  double* Uloc = calloc( N , sizeof(double));
  //double* tloc = calloc( N , sizeof(double));

  double beta[4]  = {0.0,0.5,0.5,1.0};
  double gamma[4] = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
// Init Initial Condition
  initU(U,h,N);
//for(int j = 0; j<N ;j++) printf("%d %f\n",j,U[j]);
//for(int j = 0; j<N ;j++) printf("U: %f dudx: %f\n",U[j],dudx[j]);
// Time integrating
  for(int tn = 0; tn < Ntime ; tn++){ // Time Loop
        cpyVector(U,Us,N);
        // Solving The Matrix
        for(int k = 0; k < 4 ; k++){ // RK4 Loop
            sum_Csum_Vector(Uloc,Us,du,N,beta[k]);
            //tloc[m] = t[m]  + gamma[k]*dt;
            solverFDES3(Uloc,du,c,h,dt,N);
            sum_Csum_Vector(U,U,du,N,gamma[k]);
            //t[m] = t[m] + gamma[k]*dt;
        }
  }
  //printf("ITERATION \n");
  for(int j = 0; j<N ;j++) printf("%d %f \n",j,U[j]);

free(U);
//free(t);
//free(dudx);
free(du);
free(Us);
free(Uloc);
//free(tloc);
return 1;






}
