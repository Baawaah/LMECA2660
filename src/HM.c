
#include "HM.h"

/*
* LMECA2660 Convection-Diffusion Equation
* by Thanh-Son TRAN 8116-12-00
*
*/

int main(int argc,char* argv[]){
// Init Value

/* Parameter
 *
 */
  // As for ct/L = .25 .5.75 1;
  double CFL     =   1.0 ;
  double ct_L    =   1.0;
  double Re_sig  =  40.0 ;
  int N = 128;
  int diffFlag = 0; // 0 - Convection 1 - Convection & Diffusion
  int solverType = 1;
  /* Solver Type ------
   1 - E2
   2 - S3
   3 - E4
   4 - I4
   5 - I6
  -------------------*/



  double c = 1.0;
  double L = 1.0;
  double sigma = L/32.0;
  double h = L/N;
  double dt = CFL*h /c;
  double nu;
  if(diffFlag == 1){
    nu = c*sigma/Re_sig;
  }else{
    nu = 0;
  }


  double t = 0;
  int Ntime = L*ct_L/(c*dt);
  //double tloc = 0;
  double* U    = calloc(N,sizeof(double));
  double* du   = calloc(N,sizeof(double));
  double* d2u  = calloc(N,sizeof(double));

  double* Us   = calloc(N , sizeof(double));
  double* Uloc = calloc(N , sizeof(double));
  double* Uex  = calloc(N , sizeof(double));
  double* Udif = calloc(N , sizeof(double));

  double* Qnh  = calloc(Ntime,sizeof(double));
  double* Enh  = calloc(Ntime,sizeof(double));
  double* Rnh  = calloc(Ntime,sizeof(double));


// ---Code Benchmarking-------
  struct timespec start, finish;
  double elapsed;
  clock_gettime(CLOCK_MONOTONIC, &start);
// ---------------------------

  double beta[4]  = {0.0,0.5,0.5,1.0};
  double gamma[4] = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
// Init Initial Condition
  initU(U,sigma,h,N);
  initExactU(Uex,sigma,h,c,t,L,nu,N);
// Time integrating
  for(int tn = 0; tn < Ntime ; tn++){ // Time Loop
        cpyVector(U,Us,N);
        // Solving The Matrix
        for(int k = 0; k < 4 ; k++){ // RK4 Loop
            sum_Csum_Vector(Uloc,Us,du,N,beta[k]);
            //tloc = t  + gamma[k]*dt;
            if(solverType == 1) solverFDCEE2(Uloc,du,c,h,dt,N);
            if(solverType == 2) solverFDCES3(Uloc,du,c,h,dt,N);
            if(solverType == 3) solverFDCEE4(Uloc,du,c,h,dt,N);
            if(solverType == 4) solverFDCEI4(Uloc,du,c,h,dt,N);
            if(solverType == 5) solverFDCEI6(Uloc,du,c,h,dt,N);
            if(nu != 0){
              solverFDDEE2(Uloc,d2u,nu,h,dt,N);
              sumVector(du,du,d2u,N);
              // Comment: Oui, je sais, ce n'est pas la meilleur implémentation
              //          mais, c'était pour mieux comprendre ce que je faisais.
            }
            sum_Csum_Vector(U,U,du,N,gamma[k]);
            t = t + gamma[k]*dt;
        }
        initExactU(Uex,sigma,h,c,L,t,nu,N);
        Qnh[tn] = h*sumArray(U,N);
        Enh[tn] = 1.0/2.0 *h* sumArraySquare(U,N);
        diffVector(U,Uex,N,Udif);
        Rnh[tn] = sqrt(h*sumArraySquare(Udif,N));
   }
  //printf("ITERATION \n");

//
// ---Code Benchmarking-------
  clock_gettime(CLOCK_MONOTONIC, &finish);
  elapsed = (finish.tv_sec - start.tv_sec);
  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
// ---------------------------


printf("Convection-Diffusion Simulation Code\n");
printf("by S. TRAN\n");
printf("N: %d CFL: %f ct/L: %f Re: %f\n",N,CFL,ct_L,Re_sig);
printf("Time Elapsed: %f s\n",elapsed);
printf("Solver Type: ");
if(solverType == 1) printf("Solver E2\n");
if(solverType == 2) printf("Solver S3\n");
if(solverType == 3) printf("Solver E4\n");
if(solverType == 4) printf("Solver I4\n");
if(solverType == 5) printf("Solver I6\n");

FILE* file1 = fopen("data/data1.txt","w");
if(file1 == NULL){ fprintf(stderr,"File error\n"); exit(1);}
for(int j = 0; j<N; j++) fprintf(file1,"%d %f %f\n",j,U[j],Uex[j]);
fclose(file1);
FILE* file2 = fopen("data/data2.txt","w");
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
