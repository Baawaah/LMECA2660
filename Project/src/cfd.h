/*
 * CFD.h
 * CFD General Structure and methods titles
 */
 #include <stdlib.h>

   //Domain Physical Propiety
   struct _problem {
   double CFL;
   double L;
   double H;
   double Ls;
   double Hs;
   double h;
   double dt;
   //Domain Numerical Size
   int  Nx;
   int  Ny;
   int *imax_map;
   //Domain Data of size N_x * N_y
   double **omega;
   double **psi;
   double **u;
   double **v;
 };

/*
 * @post allocate memory to vector omega,psi,u and v
 */
struct _problem* init_problem();
struct _problem* init_problem_physical        (struct _problem* Problem, double CFL, double L, double H, double Ls, double Hs, double h, double dt);
struct _problem* init_problem_numerical       (struct _problem* Problem);
struct _problem* init_problem_vector_domain   (struct _problem* Problem);
void             free_problem_vector_domain   (struct _problem* Problem);
