
#include "HM.h"

/*
* LMECA2660 Convection-Diffusion Equation
* by Thanh-Son TRAN 8116-12-00
*
*/

int main(int argc,char* argv[]){
/*  ########################################
    Parameter
   ########################################*/
  double CFL     =   1.0 ;
  double ct_L    =   1.00; int tauFile = 100; // tauFile for file writing index
  double Re_sig  =  40.0 ;
  int N = 256;
  int diffFlag = 0; // [0] - Convection  | [1] - Convection & Diffusion
  int solverType = 5;
  /* Solver Type ------
   1 - E2
   2 - S3 (Décentré O3)
   3 - E4
   4 - I4
   5 - I6
  -------------------*/
  /*  ########################################

     ########################################*/

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
              sumVector(d2u,du,du,N);
              // Comment: Oui, je sais, ce n'est pas la meilleur implémentation
              //          mais, c'était pour mieux comprendre ce que je faisais.
              //          De même que c'est aussi pour me facilité la vie pour le
              //          code vu qu'on compare plusieurs méthodes entre eux.
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

//
// ---Code Benchmarking-------
  clock_gettime(CLOCK_MONOTONIC, &finish);
  elapsed = (finish.tv_sec - start.tv_sec);
  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
// ---------------------------


printf("Convection-Diffusion Simulation Code\n");
printf("by S. TRAN\n");
printf("N: %d CFL: %f ct/L: %f Re: %f nu: %f sigma:%f \n",N,CFL,ct_L,Re_sig,nu,sigma);
printf("Time Elapsed: %f s\n",elapsed);
// Juste pour donner une idée qualitative sur la rapidité du code.
printf("Solver Type: ");
char solverName[3];
if(solverType == 1) sprintf(solverName,"E2");
if(solverType == 2) sprintf(solverName,"S3");
if(solverType == 3) sprintf(solverName,"E4");
if(solverType == 4) sprintf(solverName,"I4");
if(solverType == 5) sprintf(solverName,"I6");
printf("Solver %s\n",solverName);
char buffA[50];
char buffB[50];
if(diffFlag == 1){
  sprintf(buffA,"data/%s_%d_%d_diff_dataA.txt",solverName,N,tauFile);
  sprintf(buffB,"data/%s_%d_%d_diff_dataB.txt",solverName,N,tauFile);
}else{
  sprintf(buffA,"data/%s_%d_%d_dataA.txt",solverName,N,tauFile);
  sprintf(buffB,"data/%s_%d_%d_dataB.txt",solverName,N,tauFile);
}
FILE* file1 = fopen(buffA,"w");
if(file1 == NULL){ fprintf(stderr,"File error\n"); exit(1);}
for(int j = 0; j<N; j++) fprintf(file1,"%d %f %f\n",j,U[j],Uex[j]);
fclose(file1);
FILE* file2 = fopen(buffB,"w");
if(file2 == NULL){ fprintf(stderr,"File error\n"); exit(1);}
for(int j = 0; j<Ntime; j++) fprintf(file2,"%d %f %f %f\n",j,Qnh[j],Enh[j],Rnh[j]);
fclose(file2);
free(U);
free(Uex);
free(Udif);
//free(t);
free(du);
free(d2u);
free(Us);
free(Uloc);
//free(tloc);
free(Qnh);
free(Enh);
free(Rnh);
return 1;






}
