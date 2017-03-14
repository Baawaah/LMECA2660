#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/sysinfo.h>
#include <stddef.h>
#include <unistd.h>

void initU(double* U_value,double sigma,double h,int N){
  double Q     = 1.0;
  for(int i = 0; i < N; i++){
    U_value[i] = Q/sqrt(M_PI*sigma*sigma) * exp(-((i-N/2)*(i-N/2)*h*h)/(sigma*sigma));
  }
  //U_value[N/2] = 10.0;
  return;
}
double negaModulo(double A,double mod){
  if(A >= 0.0){ return fmod(A,mod);}
  else{ return -(fmod(-A,mod));}

}
void initExactU(double* U_value,double sigma,double h,double c,double L,double t,int N){
  double Q = 1.0;
  for(int i = 0; i < N; i++){
    double xtilde = ( i-N/2.0)*h - c*t + 1.0/2.0 * L;
    double xhat   = xtilde - floor( xtilde/L )*L;
    U_value[i] = Q/sqrt(M_PI*sigma*sigma) * exp( -( pow( xhat - L/2 ,2))/(sigma*sigma));
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
void sumVector(double* A,double* B,int size,double inc){
    for(int i = 0; i< size; i++) B[i] = A[i]+inc;
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

void solverFDEE2(double* U, double* du,double c,double h,double dt, int N){
  for(int m = 0; m < N; m++) du[m] =  dt*(-c)*(U[(m+1)%N] - U[(m+N-1)%N]) / (2.0*h);
  return;
}
void solverFDEE4(double* U, double* du,double c,double h,double dt, int N){
  for(int m = 0; m < N; m++) du[m] =  dt*(-c)*(U[(m+2)%N]+8.0*U[(m+1)%N]-8.0*U[(m-1+N)%N]+U[(m-2+N)%N])/(12.0*h);
  return;
}
void solverFDES3(double* U, double* du,double c,double h,double dt, int N){
  for(int m = 0; m < N; m++) du[m] = dt*(-c)*(2.0*U[(m+1)%N]+3.0*U[(m)%N]-6.0*U[(m-1+N)%N]+U[(m-2+N)%N])/(6.0*h);
  return;
}
