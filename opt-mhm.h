/* Global variables and function prototypes for opt-mhm
 *
 * $Revision: 2.0 $
 * $Author: tbereau $
 * $Date: 2010/08/27 $
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#ifdef OPENMP
#include <omp.h>
#endif

/* Global variables */
extern int NUM_COORD1;
extern int NUM_COORD2;
extern double TOL_FERMI;
extern double TOL_ITER;
extern int N_SIMS;
extern double UPDATE_COEFF;
extern double COORD1_MIN;
extern double COORD1_MAX;
extern double COORD1_WIDTH;
extern double COORD2_MIN;
extern double COORD2_MAX;
extern double COORD2_WIDTH;
extern int MAXFERMI;
extern double TEMP_PROB;
extern double TMIN;
extern double TMAX;
extern double TSTEP;
extern double ESTEP;
extern int BSTRAP;
extern int PARTIAL_Q_ID;
extern double PARTIAL_Q_MIN;
extern double PARTIAL_Q_MAX;
extern int PARTIAL_Q_NUM;
extern int COORD2_FLAG;
extern int PARTIAL;
extern int EBINS;
extern double EMIN;
extern double EMAX;
extern char PROBMAX_LETTER;

extern double EPS;
extern double PI;


/* GLOBAL ARRAYS */
extern char *BETA_FILE;
extern char *OUTPUT_FILE;
extern char *TEMP_AVERAGE;
extern char *TEMP_AVERAGE2;
extern char *PARTIAL_FILE;
extern char *MICRO_FILE;
extern char *F_FILE;
extern char *B_FILE;
extern double *BETAS;
extern double *K_SPRING;
extern double *X0_SPRING;
extern int *START_SAMPLING;
extern int *AUTOCOR_FACTOR;
extern int *NORM_HIST;
extern int *HIST_SIZES;
extern double **HIST;
extern double **COORD1;
extern double **COORD2;
extern double *FENERGIES;
extern double *FENERGIES_TEMP; 
extern double **PROB1;
extern double **PROB2;
extern double **B_HIST;
extern double **B_ENTROPY;
extern double **B_ERROR;
extern double **B_PROB;


#define LINESIZE 100

void readinvtemp(char *file, int umbrella_flag);
void readfile(char *file,int sim, int set_hist_boundaries);
int readfenergies(void);
void writefenergies(void);
void sighandler(int);
void self_iterative(int umbrella_flag);
void optimizedf(int umbrella_flag, int bootstrapping);
void calc_prob(void);
void calc_prob_umbrella(void);
void temp_averages(void);
void microcanonical(void);
void density_of_states(double E, double *g_E);
void b_densityofstates(double E, double *g_E);
void b_resample(void);
void b_selfiterative(int umbrella_flag);
void b_entropy(int index);
void b_umbrella(int index);
void b_coord1(int index);
void b_error(int microcanonical_flag);
void write_b_file(int microcanonical_flag);
double halfinterval(int q, int i, int left, int umbrella_flag, int bootstrapping);
double fermi(int q, int i, double x, int left, int umbrella_flag, int bootstrapping);
double init_fermi(int i, int inc_left, int umbrella_flag, int bootstrapping);






