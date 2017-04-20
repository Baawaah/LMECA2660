/*
 * CFD.h
 * CFD General Structure and methods titles
 */
 #include <stdlib.h>
 #include <math.h>

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
   int    Nx;
   int    Ny;
   int    NLs;
   int    NHs;
   int    Ntime;
   double phi;
   int *imax_map;
   //Domain Data of size N_x * N_y
   double **omega;
   double **psi;
   double **psi_s;
   double **u;
   double **v;
 };

void             print_problem_data           (struct _problem* Problem);
/* ###################################
 *  CFD Main
 * ###################################
 */
struct _problem* init_problem();
void             init_problem_physical        (struct _problem* Problem, double CFL, double L, double H, double Ls, double Hs, double h, double dt, double tmax);
void             init_problem_numerical       (struct _problem* Problem, double phi);
void             init_problem_map             (struct _problem* Problem);
void             init_problem_vector_domain   (struct _problem* Problem);
void             free_problem_vector_domain   (struct _problem* Problem);
/* ###################################
 *  CFD Numerical
 * ###################################
 */
void             boundary_psi_update          (struct _problem* Problem, double (*Q)(double) );
void             boundary_omega_update        (struct _problem* Problem);
/* ###################################
 *  CFD Test
 * ###################################
 */
double           test_Qfunc_const             (double t);
void             test_omega_domainFill        (struct _problem* Problem);
void             test_psi_boundaryFill        (struct _problem* Problem, double (*Q)(double));
void             test_omega_boundaryFill      (struct _problem* Problem);
