
#include "thomas.h"
#include "HM.h"

int main(int argc,char* argv[]){
// Init Value
  double c = 1.0;
  int N = 128;
  double CFL = 1.0;
  // As for ct/L = .25 .5.75 1;
  double ct_L = 1.0;
  double L = 1.0;
  double sigma = L/32.0;
  double h = L/N;
  double dt = CFL*h /c;
  double t = 0;
  int Ntime = L*ct_L/(c*dt);
  //double tloc = 0;
  double* U   = calloc(N    ,sizeof(double));
  double* du  = calloc(N ,sizeof(double));

  double* Us   = calloc(N , sizeof(double));
  double* Uloc = calloc(N , sizeof(double));
  double* Uex  = calloc(N , sizeof(double));
  double* Udif = calloc(N , sizeof(double));

  double* Qnh = calloc(Ntime,sizeof(double));
  double* Enh = calloc(Ntime,sizeof(double));
  double* Rnh = calloc(Ntime,sizeof(double));


  double beta[4]  = {0.0,0.5,0.5,1.0};
  double gamma[4] = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
// Init Initial Condition
  initU(U,sigma,h,N);
  initExactU(Uex,sigma,h,c,t,L,N);
//for(int j = 0; j<N ;j++) printf("%d %f\n",j,U[j]);
//for(int j = 0; j<N ;j++) printf("U: %f dudx: %f\n",U[j],dudx[j]);
// Time integrating
  for(int tn = 0; tn < Ntime ; tn++){ // Time Loop
        cpyVector(U,Us,N);
        // Solving The Matrix
        for(int k = 0; k < 4 ; k++){ // RK4 Loop
            sum_Csum_Vector(Uloc,Us,du,N,beta[k]);
            //tloc = t  + gamma[k]*dt;
            solverFDEE4(Uloc,du,c,h,dt,N);
            sum_Csum_Vector(U,U,du,N,gamma[k]);
            t = t + gamma[k]*dt;
        }
        initExactU(Uex,sigma,h,c,L,t,N);
        Qnh[tn] = h*sumArray(U,N);
        Enh[tn] = 1.0/2.0 *h* sumArraySquare(U,N);
        diffVector(U,Uex,N,Udif);
        Rnh[tn] = sqrt(h*sumArraySquare(Udif,N));
   }
  //printf("ITERATION \n");
  printf("TIME: %f N: %d NTime: %d h: %f dt: %f \n",t,N,Ntime,h,dt);

FILE* file1 = fopen("data1.txt","w");
if(file1 == NULL){ fprintf(stderr,"File error\n"); exit(1);}
for(int j = 0; j<N; j++) fprintf(file1,"%d %f %f\n",j,U[j],Uex[j]);
fclose(file1);
FILE* file2 = fopen("data2.txt","w");
if(file2 == NULL){ fprintf(stderr,"File error\n"); exit(1);}
for(int j = 0; j<Ntime; j++) fprintf(file2,"%d %f %f %f\n",j,Qnh[j],Enh[j],Rnh[j]);
fclose(file2);
free(U);
//free(t);
//free(dudx);
free(du);
free(Us);
free(Uloc);
//free(tloc);
free(Qnh);
free(Enh);
free(Rnh);
return 1;






}
