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
   double t;
   double tmax;
   //Domain Numerical Size
   int  Nx;
   int  Ny;
   int  Ntime;
   int *imax_map;
   //Domain Data of size N_x * N_y
   double **omega;
   double **psi;
   double **u;
   double **v;
 };

void             print_problem_data           (struct _problem* Problem);
/* ###################################
 *  CFD Main
 * ###################################
 */
struct _problem* init_problem();
struct _problem* init_problem_physical        (struct _problem* Problem, double CFL, double L, double H, double Ls, double Hs, double h, double dt, double tmax);
struct _problem* init_problem_numerical       (struct _problem* Problem);
struct _problem* init_problem_map             (struct _problem* Problem);
struct _problem* init_problem_vector_domain   (struct _problem* Problem);
void             free_problem_vector_domain   (struct _problem* Problem);
/* ###################################
 *  CFD Numerical
 * ###################################
 */
void             boundary_omega_update(struct _problem* Problem);
