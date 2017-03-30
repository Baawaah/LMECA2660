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

double function(double t,double u){
  return u+1;
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


int main(int argc, char *argv[]){
  struct timespec start, finish;
  double elapsed;
  clock_gettime(CLOCK_MONOTONIC, &start);

  struct u_t result = RK4_classical(0,10,0,0,0.0001,function);

  clock_gettime(CLOCK_MONOTONIC, &finish);
  elapsed = (finish.tv_sec - start.tv_sec);
  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;


  double* error = malloc(sizeof(double)*result.size);
  for(int i = 0; i < result.size ; i ++){
      error[i] = fabs( result.u[i] - (exp(result.t[i])-1) )  ;
      fprintf(stdout,"%f %f %.30lf\n",result.t[i],result.u[i],error[i]);
  }
  fprintf(stderr,"Time elapsed: %f [s] \n",elapsed);


  free(error);
  free(result.t);
  free(result.u);
  return 1;
}
