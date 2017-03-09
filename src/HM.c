#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/sysinfo.h>
#include <stddef.h>
#include <unistd.h>

struct u_t {
  int size;
  double *u;
  double *t;
};

double finiteDiff_C02(double x,double h, double (*fct)(double)){
//  Finite Difference O(h^2)
return (fct(x+h)-fct(x-h))/(2*h);
}

struct u_t RK4_classical(double t0,double t1,double u_0,double t_0,double dt_0,double (*fct)(double,double)){

  double du = 0.0;
  double dt = dt_0;
  int n = 1 + (t1-t0)/dt;
  double beta[4]  = {0.0,0.5,0.5,1.0};
  double gamma[4] = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
  double* u_s = (double*) malloc(sizeof(double)*n);
  double* t_s = (double*) malloc(sizeof(double)*n);
  //printf("Status n: %d dt: %f\n ",n,dt);
  u_s[0] = u_0;
  t_s[0] = t_0;
  for(int j=1; j < n; j ++  ){
    double u_l;
    double t_l;
    double u = u_s[j-1];
    double t = t_s[j-1];
    for(int i=0; i < 4 ; i++){
      u_l = u_s[j-1] + beta[i]*du;
      t_l = t   + beta[i]*dt;
      du  = dt * fct(t_l,u_l);
      u = u + gamma[i] * du;
      t = t + gamma[i] * dt;
      //printf("Status u: %f t: %f du:%f dt: %f gamma: %f\n ",u,t,du,dt,gamma[i]);
    }
    u_s[j] = u;
    t_s[j] = t;
  }
  struct u_t result;
  result.size = n;
  result.u = u_s;
  result.t = t_s;
  return result;
}
void initU(double* U_value){
  U_value[1] = 10.0;
  return;
}




void cpyVector(double* A,double* B,int size){
    for(int i = 0;i < size; i++) B[i] = A[i];
}

int main(int argc,char* argv[]){
// Init Value
  int N = 64;
  int Ntime = 1;
  double h = 1.0;
  double dt = 0.1;
  double* U     = calloc(N    ,sizeof(double));
  double* t     = calloc(Ntime,sizeof(double));
  double* dudx  = calloc(N    ,sizeof(double));
  double* dudt  = calloc(N    ,sizeof(double));

  //int n = 1 + (t1-t0)/dt;
  double beta[4]  = {0.0,0.5,0.5,1.0};
  double gamma[4] = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
// Init Initial Condition
  initU(U);

//for(int j = 0; j<N ;j++) printf("U: %f dudx: %f\n",U[j],dudx[j]);


// Solving The Matrix
  for(int k = 0; k<N; k++){
    dudx[k] = (U[(k+1)%N] - U[(k+N-1)%N]) / (2.0*h);
  }

// Time integrating
  for(int tn = 0; tn < Ntime ; tn++){ // Time Loop
    for(int n = 0; tn < N , n++){ // Space Loop
      double* Us = calloc( N , sizeof(double));
      for(int k = 0; k < 4 ; k++){ // RK4 Loop
        Uk = Us + beta[k]*dudt; TOODOOOOOO HEEEEREEEE
        tp = t  + gamma[k]*dt;
        dudt[n] = dt*FX TOO DOOOO;
        U = U + gamma[k]*dudt;
        t = t + gamma[k]*dt;
      }
    }
  }


for(int j = 0; j<N ;j++) printf("U: %f dudx: %f \n",U[j],dudx[j]);
printf("%d %d\n",4%4,(4+1)%4);
free(U);
free(dudx);
return 1;






}
