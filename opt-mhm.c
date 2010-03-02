/*
 * opt-mhm -- Optimized Multiple Histogram Method
 * Computes thermodynamic properties and continuous approximations
 * to different observables.
 *
 * Version: 1.7
 * Author: Tristan Bereau 
 * Date: 10/2009
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "opt-mhm.h"

#ifdef OPENMP
#include <omp.h>
#endif

// Useful macros
#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif


double COORD1_MIN, COORD1_MAX, COORD1_WIDTH;
double COORD2_MIN, COORD2_MAX, COORD2_WIDTH;
int NUM_COORD1, NUM_COORD2, N_SIMS;
double TEMP_PROB, TMIN, TMAX;
char *META_FILE;
double *BETAS, *FENERGIES, *FENERGIES_TEMP;
int *NORM_HIST, *HIST_SIZES, *START_SAMPLING, *AUTOCOR_FACTOR;
double **HIST, **COORD1, **COORD2, **PROB1, **PROB2; 
int COORD1_FLAG, COORD2_FLAG;
int PARTIAL;
int PARTIAL_Q_ID, PARTIAL_Q_NUM;
double PARTIAL_Q_MIN, PARTIAL_Q_MAX;


double EPS          =    1e-15;
double PI           = 3.14159265358979323846;
double UPDATE_COEFF = 0.5;      /* Update coefficient of the SINH algorithm. A higher value is faster,
				 * but also less stable... */
int    MAXFERMI     = (int)1e5; // Maximum number of trials before claiming the function has no solution.
double TOL_FERMI    =    1e-12; // Tolerance when solving fermi equation.
double TOL_ITER     =     1e-9; // tolerance when converging free energies.
double TSTEP        =     0.01; // Temperature step between WHAM averages
double ESTEP        =       1.; // Energy step for microcanonical analysis
char *PARTIAL_FILE  = "partial.dat";
char *OUTPUT_FILE   = "f.profile.dat";
char *TEMP_AVERAGE  = "avg_1.dat";
char *TEMP_AVERAGE2 = "avg_2.dat";
char *MICRO_FILE    = "micro.dat";
char *F_FILE        = "free_energies.dat";



// Init message
char *init = "\n\
        Optimized Multiple Histogram Method\n\
          Tristan Bereau (bereau@cmu.edu)\n\
                     2008-2009\n\
\n\
";



// Command line arguments
char *COMMAND_LINE = "options:\n\
  -d                        apply DI method *after* SINH (Bereau and Swendsen, J Comp Phys, 2009)\n\
  -do                       apply DI method ONLY (default is SINH only)\n\
  -e  energy_step           energy step for micro- and canonical analyses (default 1.)\n\
  -f  output                output file name for the free energy profile (default 'f.profile.dat')\n\
  -f1 output                output file name for WHAM temperature averages (default 'avg_1.dat')\n\
  -f2 output                output file name for canonical temperature averages (default 'avg_2.dat')\n\
  -ff file                  input/output file name for free energies (default 'free_energies.dat')\n\
  -hr temp                  analyse Hamiltonian Replica Exchange data at temperature 'temp'\n\
  -m                        turn on microcanonical analysis\n\
  -mf file                  output file for microcanonical analysis (default 'micro.dat')\n\
  -p  qid qmin qmax         calculation of one order parameter with respect to the other\n\
                            (qid:'1' or '2') in the range [qmin:qmax].\n\
  -pf file                  save partial calculation ('-p' option) using the file name 'file'\n\
                            (default 'partial.dat')\n\
  -q1 qbins qmin qmax       add an order parameter Q1 with qbins bins in the range [qmin:qmax]\n\
  -q2 qbins qmin qmax temp  add an order parameter Q2 with qbins bins in the range [qmin:qmax]\n\
                            and calculate free energy profile at temperature temp\n\
  -tb tmin tmax             temperature boundaries for canonical calculations [tmin:tmax]\n\
                            (default given by simulation temperature extrema)\n\
  -ts tstep                 temperature step in canonical calculations (default 0.01)\n\
  -u  coefficient           update coefficient for the SINH algorithm (0;1]. More is faster, but\n\
                            also less stable (default 0.5).\n\
";


/*
 * META_FILE is of the form
 *
 * path/of/file1          temperature1          [time1 at which sampling starts]
 * path/of/file2          temperature2          [time2 at which sampling starts]
 * ...
 *
 */





int main (int argc, char * const argv[]) 
{
  time_t start_n, end_n;
  double dif_n, f_bar;
  int self_it, sinh_alg, tempbound;
  int i, micro_flag, hremd_flag;

  // By default: DI is off and SINH is on
  self_it=0; 
  sinh_alg=1;
	
  COORD1_FLAG = 0;	
  COORD2_FLAG = 0;	
  PARTIAL = 0;
  tempbound = 0;
  micro_flag = 0;
  hremd_flag = 0;
	

  // output init output
  printf("%s\n",init);
   
	
  for (i = 1; i < argc; ++i){  /* Skip argv[0] (program name). */
    /* Process optional arguments. */
    if (strcmp(argv[i], "-d") == 0)  
      self_it=1;
    else if (strcmp(argv[i], "-do") == 0){
      self_it=1;
      sinh_alg=0;
    }		
    else if (strcmp(argv[i], "-p") == 0){
      // Need 3 additional arguments
      if (argc-i<4) {
	fprintf(stderr,"Not enough arguments for option %s.\n",argv[i]);
	exit(1);									
      }		   
      PARTIAL = 1;
      PARTIAL_Q_ID  = atoi(argv[i+1]);
      if (PARTIAL_Q_ID!=1 && PARTIAL_Q_ID!=2){
	fprintf(stderr,"qid must be either 1 or 2 (%d).\n",
		PARTIAL_Q_ID);
	exit(1);
      }
      PARTIAL_Q_MIN = atof(argv[i+2]);
      PARTIAL_Q_MAX = atof(argv[i+3]);
      if (PARTIAL_Q_MIN>PARTIAL_Q_MAX){
	fprintf(stderr,"qmin can't be bigger than qmax for option %s.\n",argv[i]);
	exit(1);
      }		   
      i+=3;			
    }
    else if (strcmp(argv[i], "-pf") == 0){
      ++i;			
      PARTIAL_FILE  =      argv[i];
    }
    else if (strcmp(argv[i], "-q1") == 0){
      // Need 3 additional arguments
      if (argc-i<4) {
	fprintf(stderr,"Not enough arguments for option %s.\n",argv[i]);
	exit(1);									
      }
      NUM_COORD1   = atoi(argv[i+1]);	
      COORD1_MIN   = atof(argv[i+2]);
      COORD1_MAX   = atof(argv[i+3]);
      if (NUM_COORD1 < 2) {
	printf("The number of bins needs to be bigger than 1 (%d) for option %s.\n",NUM_COORD1, argv[i]);
	exit(1);	   
      }
      if (COORD1_MIN>COORD1_MAX){
	fprintf(stderr,"qmin can't be bigger than qmax for option %s.\n",argv[i]);
	exit(1);
      }		   
      COORD1_WIDTH = (COORD1_MAX - COORD1_MIN)/((double) NUM_COORD1);	
      COORD1_FLAG = 1;			
      i+=3;				
    }
    else if (strcmp(argv[i], "-q2") == 0){
      // Need 4 additional arguments
      if (argc-i<5) {
	fprintf(stderr,"Not enough arguments for option %s.\n",argv[i]);
	exit(1);									
      }
      NUM_COORD2 = atoi(argv[i+1]);
      COORD2_MIN = atof(argv[i+2]);
      COORD2_MAX = atof(argv[i+3]);
      TEMP_PROB  = atof(argv[i+4]);
      if (NUM_COORD2 < 2) {
	printf("The number of bins needs to be bigger than 1 (%d) for option %s.\n",NUM_COORD2, argv[i]);
	exit(1);	   
      }if (COORD2_MIN>COORD2_MAX){
	fprintf(stderr,"qmin can't be bigger than qmax for option %s.\n",argv[i]);
	exit(1);
      }		   
      COORD2_WIDTH = (COORD2_MAX - COORD2_MIN)/((double) NUM_COORD2);	
      COORD2_FLAG = 1;		
      i+=4;				
    }
    else if (strcmp(argv[i], "-f") == 0){
      if (argc-i<2) {
	fprintf(stderr,"Not enough arguments for option %s.\n",argv[i]);
	exit(1);									
      }
      OUTPUT_FILE = argv[i+1];
      ++i;			
    }
    else if (strcmp(argv[i], "-f1") == 0){
      if (argc-i<2) {
	fprintf(stderr,"Not enough arguments for option %s.\n",argv[i]);
	exit(1);									
      }
      TEMP_AVERAGE  = argv[i+1];
      ++i;			
    }
    else if (strcmp(argv[i], "-f2") == 0){
      if (argc-i<2) {
	fprintf(stderr,"Not enough arguments for option %s.\n",argv[i]);
	exit(1);									
      }
      TEMP_AVERAGE2 = argv[i+1];
      ++i;			
    }
    else if (strcmp(argv[i], "-tb") == 0){
      if (argc-i<3) {
	fprintf(stderr,"Not enough arguments for option %s.\n",argv[i]);
	exit(1);									
      }
      TMIN = atof(argv[i+1]);
      TMAX = atof(argv[i+2]);
      tempbound = 1;			
      i+=2;			
    }
    else if (strcmp(argv[i], "-ts") == 0){
      if (argc-i<2) {
	fprintf(stderr,"Not enough arguments for option %s.\n",argv[i]);
	exit(1);									
      }
      TSTEP = atof(argv[i+1]);
      ++i;			
    }
    else if (strcmp(argv[i], "-u") == 0){
      if (argc-i<2) {
	fprintf(stderr,"Not enough arguments for option %s.\n",argv[i]);
	exit(1);									
      }
      UPDATE_COEFF = atof(argv[i+1]);
      if (!(UPDATE_COEFF>0. && UPDATE_COEFF<=1.)){
	fprintf(stderr,"The update coefficient for the SINH algorithm needs to be in the interval (0;1]\n");
	exit(1);				
      }
      ++i;			
    }
    else if (strcmp(argv[i], "-e") == 0){
      if (argc-i<2) {
	fprintf(stderr,"Not enough arguments for option %s.\n",argv[i]);
	exit(1);									
      }
      ESTEP = atof(argv[i+1]);
      ++i;			
    }
    else if (strcmp(argv[i], "-m") == 0)
      micro_flag = 1;			
    else if (strcmp(argv[i], "-hr") == 0) {
      if (argc-i<2) {
	fprintf(stderr,"Not enough arguments for option %s.\n",argv[i]);
	exit(1);									
      }
      TEMP_PROB = atof(argv[i+1]);
      hremd_flag = 1;
      ++i;			
    }
    else if (strcmp(argv[i], "-mf") == 0){
      if (argc-i<2) {
	fprintf(stderr,"Not enough arguments for option %s.\n",argv[i]);
	exit(1);									
      }
      MICRO_FILE = argv[i+1];
      ++i;			
    }
    else if (strcmp(argv[i], "-ff") == 0){
      if (argc-i<2) {
	fprintf(stderr,"Not enough arguments for option %s.\n",argv[i]);
	exit(1);									
      }
      F_FILE = argv[i+1];
      ++i;			
    }
    else
      META_FILE = argv[i];			
  }

	
  if (argc < 2) {		
    fprintf(stderr,"usage: %s [options] metafile\n%s\n",argv[0],COMMAND_LINE);
    exit(1);
  }
  // At this point, check that if COORD2 is defined COORD1 is defined as well
  if (COORD2_FLAG && !COORD1_FLAG){
    fprintf(stderr,"Error. Can't use option '-q2' without '-q1'.\n");
    exit(1);
  }
	

	

  if (PARTIAL){
    // Check that the second coordinate is defined
    if (!COORD2_FLAG){
      printf("Partial free energy calculation can only be activated with two order parameters.\n");
      exit(-1);
    }		
    // Calculate the number of bins in the partial set
    if (PARTIAL_Q_ID==1)
      PARTIAL_Q_NUM = (int)((PARTIAL_Q_MAX - PARTIAL_Q_MIN)/COORD1_WIDTH + 0.5);		
    else
      PARTIAL_Q_NUM = (int)((PARTIAL_Q_MAX - PARTIAL_Q_MIN)/COORD2_WIDTH + 0.5);
    // Check integrity of PARTIAL_Q boundaries
    if (PARTIAL_Q_ID==1){
      if (PARTIAL_Q_MIN < COORD1_MIN && PARTIAL_Q_MAX > COORD1_MAX) {
	printf("Partial coordinate: interval [%f,%f] out of initial range [%f,%f].\n",
	       PARTIAL_Q_MIN,PARTIAL_Q_MAX,COORD1_MIN,COORD1_MAX);
	exit(-1);				
      } else {
	if (PARTIAL_Q_MIN < COORD2_MIN && PARTIAL_Q_MAX > COORD2_MAX) {
	  printf("Partial coordinate: interval [%f,%f] out of initial range [%f,%f].\n",
		 PARTIAL_Q_MIN,PARTIAL_Q_MAX,COORD2_MIN,COORD2_MAX);
	  exit(-1);
	}
      }			
    }

  }

  /*
   * Load inverse temperatures and histogram files.
   * Stores inverse temperatures in *BETAS and
   * load all histogram files.
   */
  readinvtemp(META_FILE);
  
  /* If we're analyzing HREMD data, values of BETAS are actually coupling
   * parameters.  Subtract one from the coupling to get the contribution of
   * the altered force field only.  Make sure the energy associated with it is
   * in fact the contribution of the altered force field.
   */
  if (hremd_flag) {
    for (i=0; i<N_SIMS; ++i)
      BETAS[i] = TEMP_PROB * (BETAS[i]-1.);    
  }

  // Set Tmin and Tmax according to the lowest and highest simulations
  // .. only if the user hasn't set them explicitly
  if (!tempbound){
    TMIN=1e30;
    TMAX=-1e30;		
    for (i=0; i<N_SIMS; ++i){
      TMIN = min(TMIN,1./BETAS[i]);
      TMAX = max(TMAX,1./BETAS[i]);
    }
  }
	
  printf("Temperature averages will be calculated in the range [%2.2f,%2.2f]\n",TMIN,TMAX);
	

  // Initialize free energy arrays
  FENERGIES = calloc(N_SIMS * sizeof *FENERGIES, sizeof *FENERGIES);
  FENERGIES_TEMP = calloc(N_SIMS * sizeof *FENERGIES_TEMP, sizeof *FENERGIES_TEMP);

  // Attempt to read free energies from file
  if (readfenergies() == 0) {
    // In this case, reading the file failed and we reallocate arrays
    free(FENERGIES);
    free(FENERGIES_TEMP);		
    FENERGIES = calloc(N_SIMS * sizeof *FENERGIES, sizeof *FENERGIES);
    FENERGIES_TEMP = calloc(N_SIMS * sizeof *FENERGIES_TEMP, sizeof *FENERGIES_TEMP);
  }
	

	
  printf("\nTarget convergence: %e\n",TOL_ITER);	
	
  time (&start_n);   
  // *** Optimized MHM algorithm - SINH ***
  if (sinh_alg)
    optimizedf();
  // Self-iterative method (DI)
  if (self_it)
    self_iterative();	

  time (&end_n);
	
  dif_n=difftime(end_n,start_n);
  printf("time of computation: %d seconds.\n\n",(int)dif_n);

  // Write free energies to file
  writefenergies();
	

  // Calculate WHAM probabilities if at least one order parameter has been set
  // Note that calc_prob() can handle 1 or 2 order parameters.
  if (COORD1_FLAG)
    calc_prob();


  // Calculate averages with respect to temperature--unless we're analyzing HREMD data
  if (!hremd_flag)
    temp_averages();	


  // Shift all free energies such that <f> = 0.
  f_bar = 0.;
  for (i=0; i<N_SIMS; ++i)
    f_bar += FENERGIES[i];		
  f_bar /= N_SIMS;
  for (i=0; i<N_SIMS; ++i)
    FENERGIES[i] -= f_bar;
	
	
	
  // microcanonical analysis (optional)
  if (micro_flag)
    microcanonical();
  

	
  // Free stuff
  free(BETAS);
  free(FENERGIES);	
  free(FENERGIES_TEMP);
  free(NORM_HIST);
  if (COORD1_FLAG){
    for (i=0;i<N_SIMS;++i)
      free(COORD1[i]);
    free(COORD1);
  }
  if (COORD2_FLAG){
    for (i=0;i<N_SIMS;++i)
      free(COORD2[i]);
    free(COORD2);
  }
	
  for (i=0;i<N_SIMS;++i)
    free(HIST[i]);
  free(HIST);	
	

  printf("Operation successful.\n");   

  return 0;
}



// Read temperatures and convert to 1/T in BETAS array
void readinvtemp(char *file)
{
  FILE *ffile;
  char *line;
  char path_i[LINESIZE];
  int i=0, vals, sampling_i, autocor_i;
  double temp_i;

  // Create line.
  line = (char *) malloc(sizeof(char) * LINESIZE);
  if (!line)
    {
      printf("couldn't allocate space for line\n");
      exit(-1);
    }


  // Read a first time to obtain the number of simulations
  // and the histogram boundaries.
  ffile = fopen(file, "r");	
  if (!ffile)
    {
      free(line);
      printf("cannot open file %s.\n",file);	   
      exit(-1);		
    }

  HIST_SIZES = NULL;	
  BETAS = NULL;
  START_SAMPLING = NULL;   
  AUTOCOR_FACTOR = NULL;   
	
  line = fgets(line,LINESIZE,ffile);	
  while (line != NULL) {
    if (line[0] !=  '#') {
      vals = sscanf(line,"%s %lf %d %d", path_i, &temp_i, &sampling_i,&autocor_i);
      if (!(vals > 1 && vals < 5)) {
	printf("failure reading %s : can't read (path, temp, [sampling_start, autocor_factor])\n", file);
	exit(-1);				
      }
      if (vals < 3)		   	
	sampling_i = 0;
      if (vals < 4)
	autocor_i = 0;				
			
			
      // Grow arrays dynamically
      HIST_SIZES     = (int*) realloc (HIST_SIZES, (i+1)*sizeof(int));			   
      BETAS          = (double*) realloc (BETAS, (i+1)*sizeof(double));
      START_SAMPLING = (int*) realloc (START_SAMPLING, (i+1)*sizeof(int));
      AUTOCOR_FACTOR = (int*) realloc (AUTOCOR_FACTOR, (i+1)*sizeof(int));			
      // Initialize HIST_SIZES
      HIST_SIZES[i] = 0;				
      // Store \beta=1/T (in units of k_B*T_r).
      // Note that the temperature is read and converted to 1/T
      BETAS[i] = 1./temp_i;
      // sampling starts at sampling_i
      START_SAMPLING[i] = sampling_i;
      // autocorrelation factor for simulation i
      AUTOCOR_FACTOR[i] = autocor_i;		   
			
      // Determine histogram boundaries
      readfile(path_i,i,1);				
			
      ++i;
			
      line = fgets(line,LINESIZE,ffile);	 
    } else {
      line = fgets(line,LINESIZE,ffile);
    }		
  }
  fclose(ffile);
  free(line);

  // From the number of temperatures entered, determine N_SIMS.
  N_SIMS = i;
  // Restart counter
  i = 0;
  // number of hits for each histogram
  NORM_HIST = calloc (N_SIMS * sizeof *NORM_HIST, sizeof *NORM_HIST);
  // Now read again to go through the list of histograms.
  ffile = fopen(file, "r");

  // Create line.
  line = (char *) malloc(sizeof(char) * LINESIZE);
  if (!line)
    {
      printf("couldn't allocate space for line\n");
      exit(-1);
    }

  if (!ffile)
    {
      free(line);
      //fclose(ffile);
      printf("cannot open file %s.\n",file);
      exit(-1);
    }
  line = fgets(line,LINESIZE,ffile);
  while (line != NULL) {
    if (line[0] !=  '#') {
      vals = sscanf(line,"%s %lf", path_i, &temp_i);
      if (vals > 2) {
	printf("failure reading %s : can't read (path temp)\n", file);
	exit(-1);
      }
      if (vals == 2) {
	// Initialize NORM_HIST
	NORM_HIST[i] = 0;				

	// Load histograms
	readfile(path_i,i,0);
	++i;
      }
      line = fgets(line,LINESIZE,ffile);
    } else {
      line = fgets(line,LINESIZE,ffile);
    }
  }
  fclose(ffile);


  free(line);



}



void readfile(char *histo_file, int sim, int set_hist_boundaries)
{
  char *line;
  FILE *file;
  int vals, j, counter, k_autocor_counter;
  double time, value, value2, energy;
	
  counter=0;
	
  if (set_hist_boundaries==0){
    if (sim==0) {
      // Energy histogram
      HIST   = calloc (N_SIMS * sizeof **HIST, sizeof **HIST);			
      if (HIST==NULL){
	printf("Pointer was not initialized properly at line %d.\n",__LINE__);
	exit(-1);
      }
      // 1st order parameter
      if (COORD1_FLAG){
	COORD1 = calloc (N_SIMS * sizeof **COORD1, sizeof **COORD1);
	if (COORD1==NULL){
	  printf("Pointer was not initialized properly at line %d.\n",__LINE__);
	  exit(-1);
	}
      }
      // 2nd order parameter				
      if (COORD2_FLAG){				
	COORD2 = calloc (N_SIMS * sizeof **COORD2, sizeof **COORD2);
	if (COORD2==NULL){
	  printf("Pointer was not initialized properly at line %d.\n",__LINE__);
	  exit(-1);
	}
      }		   
      for (j = 0; j < N_SIMS; ++j) {
	HIST[j]   = calloc(HIST_SIZES[j] * sizeof *HIST, sizeof *HIST);
	if (COORD1_FLAG)
	  COORD1[j] = calloc(HIST_SIZES[j] * sizeof *COORD1, sizeof *COORD1);
	if (COORD2_FLAG)
	  COORD2[j] = calloc (HIST_SIZES[j] * sizeof *COORD2, sizeof *COORD2);
      }
    }
  } 
	
  line = malloc(sizeof(char) * LINESIZE);
  if (!line)
    {
      fprintf(stderr,"couldn't allocate space for line\n");
      exit(1);
    }
	
  file = fopen(histo_file, "r");
  if (!file)
    {
      free(line);
      //fclose(file);
      printf("cannot open file %s.\n",histo_file);	   
      exit(-1);
    }
	
  line = fgets(line,LINESIZE,file);
	
  // Initialize the autocorrelation counter to 0
  k_autocor_counter = 0;
	
  while (line != NULL)
    {
      if (line[0] != '#')
	{
	  /* Read file of the form :
	     time energy
	     OR
	     time value1 energy
	     OR
	     time value1 value2 energy
	     where value corresponds to an order parameter
	  */
	  if (COORD2_FLAG){
	    vals=sscanf(line,"%lf %lf %lf %lf", &time, &value, &value2, &energy);
	    if (vals!=4){
	      printf("failure reading %s.\n",histo_file);
	      exit(-1);			
	    }
	  } else if (COORD1_FLAG) {
	    vals=sscanf(line,"%lf %lf %lf", &time, &value, &energy);
	    if (vals!=3){
	      printf("failure reading %s.\n",histo_file);
	      exit(-1);			
	    }
	  } else {
	    vals=sscanf(line,"%lf %lf", &time, &energy);
	    if (vals!=2){
	      printf("failure reading %s.\n",histo_file);
	      exit(-1);			
	    }
	  }
			
	  if (set_hist_boundaries==0) {
	    if (time >= START_SAMPLING[sim]){
	      if (AUTOCOR_FACTOR[sim] > 0){
		if (k_autocor_counter % AUTOCOR_FACTOR[sim]==0){
		  HIST[sim][counter]        = energy;
		  if (COORD1_FLAG)
		    COORD1[sim][counter]      = value;
		  if (COORD2_FLAG)
		    COORD2[sim][counter]  = value2;					
		  NORM_HIST[sim]           += 1;
		  ++counter;
		}
	      } else {
		HIST[sim][counter]        = energy;
		if (COORD1_FLAG)
		  COORD1[sim][counter]      = value;
		if (COORD2_FLAG)
		  COORD2[sim][counter]  = value2;					
		NORM_HIST[sim]           += 1;
		++counter;
	      }						
	    }				
	  } else {
	    if (time >= START_SAMPLING[sim]){
	      if (AUTOCOR_FACTOR[sim] > 0){						
		if (k_autocor_counter % AUTOCOR_FACTOR[sim]==0)
		  ++HIST_SIZES[sim];
	      } else
		++HIST_SIZES[sim];					
	    }				
	  }

	  // Increment the autocorelation counter
	  ++k_autocor_counter;			
	}
      line = fgets(line,LINESIZE,file);
    }

  fclose(file);
}



void optimizedf()
{
  double converg_rate, sumF;
  int iter=0, i, q, iter_q;

	
  converg_rate=1.;
  sumF=0.;	
	
  printf("\nStarting optimized MHM algorithm : \n");

  // IMPLEMENT EARLY EXIT CONDITION IF MORE NEIGHBORS (HIGHER Q) DON'T CONTRIBUTE.
  for (q=0;q<N_SIMS-1;++q) {
    printf("Calculation neighbor level q = %d\n",q);		
    converg_rate=1.;		
    iter_q = 0;		
    while (converg_rate>TOL_ITER){
      converg_rate=0.;
#ifdef OPENMP
#pragma omp parallel for reduction(+:converg_rate)
#endif
      for (i=0;i<N_SIMS-1;++i){
	if (iter==0 && q==0) {
	  FENERGIES[i] = halfinterval(q,i,1);
	} else {
	  FENERGIES_TEMP[i] = FENERGIES[i];
	  // weight the new solution in order to avoid unstability
	  FENERGIES[i] = UPDATE_COEFF*halfinterval(q,i,1) + 
	    (1-UPDATE_COEFF)*FENERGIES_TEMP[i];
	}
	converg_rate += fabs(FENERGIES[i]-FENERGIES_TEMP[i]);
      }

      printf("iteration: %d; absolute error: %e\n",iter,converg_rate);
      // Test early exit condition (more neighbors do not contribute)
      if (q > 0 && converg_rate < TOL_ITER && iter_q == 0)
	q = N_SIMS;
      ++iter;
      ++iter_q;		 
    }
		
  }
	
	
  /*
   * Transform free energy differences to free energies, and output.
   */
	
  printf("\nFinal Answers for reduced free energies :\n");
  FENERGIES[0]= 0.;	
  for (i=0;i<N_SIMS-1;++i){
    sumF += FENERGIES_TEMP[i];	   
    FENERGIES[i+1] = sumF;
    FENERGIES_TEMP[i]=FENERGIES[i];		
  }
  for (i=0;i<N_SIMS;++i)
    printf("f_%d:\t%lf\n",i+1,FENERGIES[i]);
}

int readfenergies(void)
// Attempt to read free energies from a file. Will speed up convergence for pre-existing data
{
  char *line;	
  FILE *file;
  int vals, counter;
  double value;
	
  line = malloc(sizeof(char) * LINESIZE);
  if (!line) {
    fprintf(stderr,"couldn't allocate space for line\n");
    exit(1);
  }

  file = fopen(F_FILE, "r");
  if (!file) {
    free(line);
    printf("Can't open %s.\n",F_FILE);
    return 0;		
  }

  counter = 0;	
		
  line = fgets(line,LINESIZE,file);
  while (line != NULL){
    if (line[0] != '#'){
      vals=sscanf(line,"%lf",&value);			
      if (vals!=1){
	printf("failure reading %s.\n",F_FILE);
	return 0;
      }
      FENERGIES[counter] = value;
      FENERGIES_TEMP[counter] = value;			
      ++counter;			
    }
    line = fgets(line,LINESIZE,file);		
  }
  fclose(file);

  if (counter == N_SIMS){
    printf("File %s was successfully read.\n",F_FILE);		
    return 1;
  } else {
    printf("Failed to load free energies from file %s.\n",F_FILE);
    printf("All free energies will be initialized to 0.\n");		
    return 0;
  }
	
}

void writefenergies(void)
// Write free energies to file
{
  int i;   
  FILE *file;

  printf("Saving free energies to output file %s.\n",F_FILE);
	
  file = fopen(F_FILE, "wt");
  if (file) {
    for(i = 0; i<N_SIMS; ++i) 
      fprintf(file,"%4.16f\n",FENERGIES[i]);		
  } else 
    fprintf(stderr,"Failed to output free energies to file %s.\n",
	    F_FILE);		
	
  fclose(file);
}


void self_iterative(void)
{
  double *fold_rec, *argarray, deltaF;
  double sumNum, sumDen, arg;	
  int iter, i_HE, i, j, k;

  deltaF=1.;
	
	
  fold_rec = calloc (N_SIMS * sizeof *fold_rec, sizeof *fold_rec);

  for (i = 0; i<N_SIMS; ++i)
    FENERGIES_TEMP[i] = FENERGIES[i];

  iter=0;
	

  printf("\nStarting self-iterative algorithm : \n");
  while (deltaF>TOL_ITER) {
#ifdef OPENMP
#pragma omp parallel for private(j,i_HE,sumNum,sumDen,arg,k,argarray)
#endif
    for (i = 0; i<N_SIMS; ++i){
      FENERGIES[i] = 0.;
      for (j = 0; j<N_SIMS; ++j){
	for (i_HE = 0; i_HE<HIST_SIZES[j]; ++i_HE){
	  sumDen = 0.;
	  arg    = -1e300;
	  argarray = calloc (N_SIMS * sizeof *argarray, sizeof *argarray);
	  // Determine largest argument (overflow trick)
	  for (k = 0; k<N_SIMS; ++k){
	    argarray[k] = (BETAS[i]-BETAS[k])*HIST[j][i_HE]-FENERGIES_TEMP[k];
	    // calculate max value
	    if (argarray[k]>arg)
	      arg=argarray[k];
	  }
	  // Now perform the calculation, by using the overflow trick
	  for (k = 0; k<N_SIMS; ++k)
	    sumDen += NORM_HIST[k]*exp(argarray[k]-arg);
	  sumNum = exp(-arg);
	  FENERGIES[i] += sumNum/sumDen;
	  free(argarray);
	}
      }
      FENERGIES[i]      = log(FENERGIES[i]);
      fold_rec[i]       = FENERGIES_TEMP[i];
      FENERGIES_TEMP[i] = FENERGIES[i];
    }
    deltaF = 0.;
    for (i = 0; i<N_SIMS; ++i )
      deltaF += fabs(FENERGIES[i]-fold_rec[i]);
    ++iter;
    printf("%d \t delta : %e\n",iter,deltaF);

  }

  printf("\n");
  for (i = 0; i<N_SIMS; ++i)
    printf("%d \t %f\n",i,FENERGIES[i]);	

  free(fold_rec);
	
}


void calc_prob(void) 
{
  int m, n, i, j, i_HE, k, t_count;
  FILE *file;   
  double *argarray;	
  double sumNum, sumDen, arg, bin1_min, bin1_max, bin2_min, bin2_max, max_prob, coord1ij, coord2ij;	
  double part_min, part_num_bar, part_min_bar, part_width_bar, part_max;	
  double temp_index, temp_min, temp_max, temp_step;	
  int temp_bins;
  // double part_width;

  temp_min=TMIN;
  temp_max=TMAX;
  temp_step=TSTEP;
	
  printf("Calculating free energies as a function of the order parameter(s).\n");
	
	
  // Calculate the dimension of the PROB1 array (dim. of Q * dim. of T)
  // Nothing deep behind the '2' in front, just trying to avoid
  // segmentation faults!
  if (!COORD2_FLAG)
    temp_bins = 2+(int)((temp_max-temp_min)/temp_step);
  else
    temp_bins = 0;
	
	
  if (COORD2_FLAG){		
    // 2D Probability
    PROB2        = calloc (NUM_COORD1 * sizeof **PROB2, sizeof **PROB2);
    for (i=0;i<NUM_COORD1;++i)
      PROB2[i] = calloc (NUM_COORD2 * sizeof PROB2, sizeof * PROB2);
		
    // If partial is on, calculate 1D Probability
    if (PARTIAL){
      PROB1    = calloc (1 * sizeof **PROB1, sizeof **PROB1);
      PROB1[0] = calloc (PARTIAL_Q_NUM * sizeof *PROB1, sizeof *PROB1);
    }
		
    // Do not loop through temperature
    temp_min  = TEMP_PROB;
    temp_max  = TEMP_PROB+1.;
    temp_step = 2.;		
  } else {
    PROB1 = calloc (temp_bins * sizeof **PROB1, sizeof **PROB1);
    for (i=0;i<temp_bins;++i)
      PROB1[i] = calloc (NUM_COORD1 * sizeof *PROB1, sizeof *PROB1);
    // This will make sure the for loops below go smoothly
    NUM_COORD2   = 1;
    COORD2_MIN   = 0;
    COORD2_MAX   = 0;
    COORD2_WIDTH = 0;
  }
	
  argarray = calloc (N_SIMS * sizeof *argarray, sizeof *argarray);

  max_prob = 0.;   

  //output free energies
  file = fopen(OUTPUT_FILE, "wt");
  if (!COORD2_FLAG)
    fprintf (file,"#Free energy profile as a function of temperature and order parameter\n");		

  t_count=0;	
  for (temp_index = temp_min; temp_index <= temp_max+1e-8; temp_index += temp_step) {
    for (m = 0; m<NUM_COORD1; ++m){
      bin1_min = COORD1_MIN +  m   *COORD1_WIDTH;
      bin1_max = COORD1_MIN + (m+1)*COORD1_WIDTH;
      for (n = 0; n<NUM_COORD2; ++n){
	bin2_min = COORD2_MIN +  n   * COORD2_WIDTH;
	bin2_max = COORD2_MIN + (n+1)* COORD2_WIDTH;
	if (COORD2_FLAG)
	  PROB2[m][n] = 0.;
	else
	  PROB1[t_count][m] = 0.;
	for (i = 0; i<N_SIMS; ++i){
	  for (i_HE = 0; i_HE<HIST_SIZES[i]; ++i_HE){
	    // Numerator - only if COORD{1,2}[i][i_HE] is/are in the correct bin
	    // If Q2 is not defined, condition is automatically satisfied (0<=0...)
	    if (COORD2_FLAG)
	      coord2ij = COORD2[i][i_HE];
	    else
	      coord2ij = 0.;					
	    if (COORD1[i][i_HE] >= bin1_min && COORD1[i][i_HE] <= bin1_max &&
		coord2ij        >= bin2_min && coord2ij        <= bin2_max ){
	      sumDen = 0.;
	      arg = -1e300;
	      for (k = 0; k<N_SIMS; ++k){
		argarray[k] = (1./temp_index-BETAS[k])*HIST[i][i_HE]-FENERGIES[k];
		// calculate max value
		if (argarray[k]>arg)
		  arg = argarray[k];
	      }
	      for (k = 0; k<N_SIMS; ++k)
		sumDen += NORM_HIST[k]*exp(argarray[k]-arg);
	      sumNum = exp(-arg);
	      if (COORD2_FLAG)
		PROB2[m][n] += sumNum /sumDen;
	      else
		PROB1[t_count][m]    += sumNum / sumDen;					
	    }
	  }
	}
	if (COORD2_FLAG)
	  if (PROB2[m][n] > max_prob)
	    max_prob = PROB2[m][n];
      }
      if (!COORD2_FLAG) 
	if (PROB1[t_count][m] > max_prob)
	  max_prob = PROB1[t_count][m];
    }
    ++t_count;		
  }

  if (file) {
    if (COORD2_FLAG){
      fprintf (file,"#Free energies as a function of the order parameter(s) at temperature %f\n",temp_min);
      for (i=0;i<NUM_COORD1;++i) {
	for (j=0;j<NUM_COORD2;++j)
	  fprintf (file,"%f\t%f\t%f\n",COORD1_MIN + (i+.5)*COORD1_WIDTH,
		   COORD2_MIN + (j+.5)*COORD2_WIDTH,
		   -TEMP_PROB*log(PROB2[i][j]/max_prob));
	fprintf (file,"\n");				
      }			
    } else {
      if (!COORD2_FLAG){
	t_count=0;			
	for (temp_index = temp_min; temp_index <= temp_max+1e-8; temp_index += temp_step) {
	  for (i=0;i<NUM_COORD1;++i) {
	    if (PROB1[t_count][i]<EPS)
	      fprintf (file,"#%f\t%f\tinf\n",temp_index,COORD1_MIN + (i+.5)*COORD1_WIDTH);
	    else
	      fprintf (file,"%f\t%f\t%f\n",temp_index,COORD1_MIN + (i+.5)*COORD1_WIDTH,-temp_index*log(PROB1[t_count][i]/max_prob));
	  }			
	  fprintf (file,"\n");
	  ++t_count;				
	}
      }
    }
  }
	
	
	
  fclose(file);	

  // Now calculate partial probability
  if (PARTIAL){
    printf("\n Probability calculation of the partial Q region:\n");		
    if (PARTIAL_Q_ID==1){
      //part_width     = COORD1_WIDTH;
      part_min       = PARTIAL_Q_MIN;
      part_min_bar   = COORD2_MIN;
      part_max       = PARTIAL_Q_MAX;			
      part_width_bar = COORD2_WIDTH;			
      part_num_bar   = NUM_COORD2;			
    } else {
      //part_width     = COORD2_WIDTH;
      part_min       = PARTIAL_Q_MIN;
      part_min_bar   = COORD1_MIN;
      part_max       = PARTIAL_Q_MAX;			
      part_width_bar = COORD1_WIDTH;			
      part_num_bar   = NUM_COORD1;			
    }
    max_prob = 0.;		
    for (m=0;m<part_num_bar;++m){
      bin1_min = part_min_bar +  m   * part_width_bar;
      bin1_max = part_min_bar + (m+1)* part_width_bar;
			
      PROB1[0][m] = 0;
      for (i=0; i<N_SIMS; ++i){
	for (i_HE = 0; i_HE<HIST_SIZES[i]; ++i_HE){
	  // Numerator only if COORD{1,2}[i][i_HE] falls into the interval
	  // bin1_min bin1_max
	  if (PARTIAL_Q_ID==1){
	    coord1ij = COORD2[i][i_HE];
	    coord2ij = COORD1[i][i_HE];
	  } else {
	    coord1ij = COORD1[i][i_HE];						
	    coord2ij = COORD2[i][i_HE];
	  }
	  if (coord1ij >= bin1_min && coord1ij <= bin1_max &&
	      coord2ij >= part_min && coord2ij <= part_max) {
	    sumDen = 0.;
	    arg = -1e300;
	    for (k = 0; k<N_SIMS; ++k){
	      argarray[k] = (1./TEMP_PROB-BETAS[k])*HIST[i][i_HE]-FENERGIES[k];
	      // calculate max value
	      if (argarray[k]>arg)
		arg = argarray[k];
	    }
	    for (k = 0; k<N_SIMS; ++k)
	      sumDen += NORM_HIST[k]*exp(argarray[k]-arg);
	    sumNum = exp(-arg);
	    PROB1[0][m]   += sumNum / sumDen;										
	  }
	}				
      }
      if (PROB1[0][m] > max_prob)
	max_prob = PROB1[0][m];
      printf ("%2.5f\t%2.5f\n",.5*(bin1_min+bin1_max),PROB1[0][m]);					
    }
    //output free energies for partial
    file = fopen(PARTIAL_FILE, "wt");
		
    if (file) {
      fprintf (file,"#Free energies as a function of the order parameter over subset of other coordinate. Temperature %f\n",TEMP_PROB);
      for (i=0;i<part_num_bar;++i) {
	fprintf (file,"%lf\t%f\n",part_min_bar + (i+.5)*part_width_bar,
		 -TEMP_PROB*log(PROB1[0][i]/max_prob));
      }
      fclose(file);
    }		
  }

  // FREE the memory of PROB1 and PROB2
  if (!PARTIAL) {
    for (i=0;i<temp_bins;++i)
      free(PROB1[i]);
    free(PROB1);
  } else {
    free(PROB1[0]);
    free(PROB1);
  }
  if (COORD2_FLAG){		
    for (i=0;i<NUM_COORD1;++i)
      free(PROB2[i]);
    free(PROB2);	
  }

  free(argarray);
	
}



void temp_averages(void)
{
  double temp, entropy, spec_heat, spec_heat_old;
  double dcoorddT, coordE;
  double f_temp, arg, beta, sumDen, E, E2, E4, tmp, coord;
  double *argarray;	
  int i_HE, i, j;	
  FILE *file, *file2;

  temp=TMIN;
  entropy=0.;
  spec_heat=0.;
  dcoorddT = 0.;
  spec_heat_old=0.;

  printf("Performing temperature averages.\n");	

  argarray = calloc (N_SIMS * sizeof *argarray, sizeof *argarray);
	
  file = fopen(TEMP_AVERAGE, "wt");   

  if (file) {
    fprintf (file,"# Averages as a function of temperature\n");
    if (COORD1_FLAG)
      fprintf (file,"# T \t E \t E^2 \t E^4 \t C_V \t S \t Order parameter \t d(Q)/dT\n");	   
    else
      fprintf (file,"# T \t E \t E^2 \t E^4 \t C_V \t S\n");	   
    while (temp<=TMAX+1e-8) {
      // energy
      E     = 0.;
      // energy square
      E2    = 0.;
      // energy to the four
      E4    = 0.;
      // order parameter
      coord = 0.;
      coordE = 0.;
			
      // Calculate free energy at temperature temp.
      f_temp = 0.;
      arg = -1e300;
      beta = 1./temp;
      for (i = 0; i <N_SIMS; ++i) {
	for (i_HE = 0; i_HE<HIST_SIZES[i]; ++i_HE){
	  sumDen = 0.;
	  for (j = 0; j<N_SIMS; ++j){
	    argarray[j] = (beta-BETAS[j])*HIST[i][i_HE]-FENERGIES[j];
	    // Calculate max value.
	    if (argarray[j]>arg)
	      arg=argarray[j];
	  }
	  for (j = 0; j<N_SIMS; ++j)
	    sumDen += NORM_HIST[j]*exp(argarray[j]-arg);
	  tmp     = exp(-arg)/sumDen;
	  f_temp +=                                        tmp;
	  E      +=                        HIST[i][i_HE] * tmp;
	  E2     +=        HIST[i][i_HE] * HIST[i][i_HE] * tmp;
	  E4     += HIST[i][i_HE] * HIST[i][i_HE] 
	    * HIST[i][i_HE] * HIST[i][i_HE] * tmp;
	  if (COORD1_FLAG) {
	    coord  +=                    COORD1[i][i_HE] * tmp;
	    coordE +=      COORD1[i][i_HE]*HIST[i][i_HE] * tmp;	 
	  }
	}
      }			
      E        /=      f_temp;
      E2       /=      f_temp;
      E4       /=      f_temp;
      if (COORD1_FLAG) {
	coord    /=      f_temp;		   
	coordE   /=      f_temp;
	// derivative of order parameter wrt temperature
	dcoorddT  = beta * beta * (coordE - coord * E);
      }
      // specific heat
      spec_heat = beta * beta * (E2 - E * E);
      // Trapezoidal rule of integration 
      entropy  += .5 * TSTEP * (spec_heat+spec_heat_old) * beta;

      if (COORD1_FLAG) 
	fprintf(file,"%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f\n",
		temp,E,E2,E4,spec_heat,entropy,coord,dcoorddT);
      else
	fprintf(file,"%f \t %f \t %f \t %f \t %f \t %f\n",
		temp,E,E2,E4,spec_heat,entropy);

      spec_heat_old = spec_heat;			
      temp += TSTEP;			
    }
  }

  fclose(file);

  if (COORD1_FLAG) {
    file2 = fopen(TEMP_AVERAGE2, "wt");
    
    
    if (file2) {
      
      /* Now check the free energy calculation by computing averages at the sampled
       * temperatures for the order parameter. This does not use free energies, \
       * it's the conventional way of calculating averages in MC/MD simulations.
       */
      for (i = 0; i<N_SIMS; ++i) {
	E = 0.;			
	sumDen = 0.;
	tmp = 0.;
	for (i_HE = 0; i_HE<HIST_SIZES[i]; ++i_HE){
	  tmp    += COORD1[i][i_HE];
	}
	fprintf (file2,"%f\t%f\n",1./BETAS[i],tmp/NORM_HIST[i]);			
      }
      
      
    }
    fclose(file2);
  }
	
}


void microcanonical(void)
{
  int i, i_HE, ebins, e_index;	
  double energy, emin, emax;
  double g_E, g_E_previous, G_E;
  double *gE, *lnGE;	
  FILE *file;

  file = fopen(MICRO_FILE, "wt");

  printf("Performing microcanonical analysis and saving to file '%s'.\n",MICRO_FILE);	


  // Determine emin and emax from the histograms.
  emin= 1e30;
  emax=-1e30;		

  for (i=0; i<N_SIMS; ++i){
    for (i_HE = 0; i_HE<HIST_SIZES[i]; ++i_HE){
      emin = min(emin,HIST[i][i_HE]);
      emax = max(emax,HIST[i][i_HE]);
    }
  }
  ebins = 2+(int)((emax - emin)/ESTEP);	
	
  // initialize the phase-space volume G_E
  G_E = 0.;
  g_E_previous = 0.;	
	
  // initialize arrays
  gE = calloc (ebins * sizeof *gE, sizeof *gE);
  lnGE = calloc (ebins * sizeof *lnGE, sizeof *lnGE);

	
  if (file) {
    fprintf (file,"# Microcanonical analysis\n");		
    fprintf (file,"# Averages as a function of energy\n");
    fprintf (file,"# E\tDoS\tEntropy\tHertz entropy\tdS/dE\tC_v\n");

    energy = emin;
    e_index = 0;
    // DoS
    printf("Calculating density of states...\n");
    while (energy<=emax+1e-8) {
      density_of_states(energy,&g_E);
		    
      gE[e_index] = g_E;
		    
      // Trapezoidal rule of integration -- integrated g(E)
      G_E += .5 * ESTEP * (g_E + g_E_previous);
		    
      // Hertz entropy log(G_E)
      lnGE[e_index]  = log(G_E);
		    
      // update variables
      energy  += ESTEP;
      ++e_index;			
      g_E_previous = g_E;
    }		

    // write to file
    printf("Writing to file...\n");
    energy = emin+ESTEP/2.;
    e_index = 0;	   
    while (energy<=emax+1e-8) {
      fprintf(file, "%f \t %e \t %e \t %e\n",energy,gE[e_index],log(gE[e_index]),lnGE[e_index]);
      // update variables
      energy  += ESTEP;
      ++e_index;			
    }

		
    // close
    fclose(file);
  } else 
    fprintf(stderr,"Warning: Can't open file %s.\n",MICRO_FILE);

}


void density_of_states(double E, double *g_E)
/* Evaluate the density of states at energy E
 * where E is in [E-ESTEP/2;E+ESTEP/2]
 * Returns for a given E interval:
 * g_E the density of states
 */
{
  int i, k, i_HE;
  double sumDen, sumNum, arg, *argarray;	
  //double dgdE, d2gdE2, d2SdE2;

  argarray = calloc (N_SIMS * sizeof *argarray, sizeof *argarray);
  sumDen = 0.;
  sumNum = 0.;
	
  *g_E = 0.;	
	
  for (i = 0; i<N_SIMS; ++i){
    for (i_HE = 0; i_HE<HIST_SIZES[i]; ++i_HE){
      if (HIST[i][i_HE] >= E-ESTEP/2. && HIST[i][i_HE] <= E+ESTEP/2.){				
	sumDen = 0.;
	arg = -1e300;
	for (k = 0; k<N_SIMS; ++k){
	  argarray[k] = -BETAS[k]*HIST[i][i_HE]-FENERGIES[k];
	  // calculate max value
	  if (argarray[k]>arg)
	    arg = argarray[k];
	}
	sumNum = exp(-arg);				
	for (k = 0; k<N_SIMS; ++k)
	  sumDen += NORM_HIST[k]*exp(argarray[k]-arg);
	*g_E += sumNum/sumDen;				
      }
    }
  }		
}


double halfinterval(int q, int i, int left)
{
  double a, b, fa, fb, d, fd;
	
  // Find initial conditions of false position method by using Bennett's equation
  a=init_fermi(i, 1);
  b=init_fermi(i, 0);


  // ***** False position Method
  do{
    fa=fermi(q, i, a, left);
    fb=fermi(q, i, b, left);
    d=(fb*a-fa*b)/(fb-fa);	   
    fd=fermi(q, i, d, left);		
    if (fa*fd>0)
      a=d;
    else
      b=d;
  } while (fabs(fd) > TOL_FERMI);
  return d;

	 
}


double fermi(int q, int i, double x, int left)
{
  double func;
  double dbel, dber, den, dbm;
  double deltafk;
  int n, i_HE, m, k;

  func=0.;
	
	
  for (n=i-q;n<=i+1+q;++n){	
    if (n>=0 && n<N_SIMS){
      for (i_HE=0; i_HE<HIST_SIZES[n]; ++i_HE){
	dbel=(BETAS[i]-BETAS[i+1])*HIST[n][i_HE];		
	dber=(BETAS[i+1]-BETAS[i])*HIST[n][i_HE];		
	den=1;
	if (left){ // coming from the left
	  den=NORM_HIST[i]+NORM_HIST[i+1]*exp(dbel-x);
	  if (q>0){
	    for (m=i-q;m<=i-1;++m){
	      if (m>=0){
		dbm=(BETAS[i]-BETAS[m])*HIST[n][i_HE];						
		deltafk=0.;
		for (k=m;k<=i-1;++k){
		  deltafk+=FENERGIES[k];
		}
		den+=NORM_HIST[m]*exp(dbm+deltafk);
	      }
	    }
	    for  (m=i+2;m<=i+1+q;++m){
	      if (m<N_SIMS){
		dbm=(BETAS[i]-BETAS[m])*HIST[n][i_HE];						
		deltafk=0.;
		for (k=i+1;k<=m-1;++k){
		  deltafk+=FENERGIES[k];
		}
		den+=NORM_HIST[m]*exp(dbm-deltafk-x);
	      }
	    }
	  }
	}
	else{ //coming from the right
	  den=NORM_HIST[i]*exp(dber+x)+NORM_HIST[i+1];
	  if (q>0){
	    for (m=i-q;m<=i-1;++m){
	      if (m>=0){
		dbm=(BETAS[i+1]-BETAS[m])*HIST[n][i_HE];						
		deltafk=0.;
		for (k=m;k<=i-1;++k){
		  deltafk+=FENERGIES[k];
		}
		den+=NORM_HIST[m]*exp(dbm+deltafk+x);
	      }
	    }
	    for  (m=i+2;m<=i+1+q;++m){
	      if (m<N_SIMS){
		dbm=(BETAS[i+1]-BETAS[m])*HIST[n][i_HE];						
		deltafk=0.;
		for (k=i+1;k<=m-1;++k){
		  deltafk+=FENERGIES[k];
		}
		den+=NORM_HIST[m]*exp(dbm-deltafk);
	      }
	    }
	  }
	}
	func+=1./den;								
      }
    }
  }
  func-=1.;
  return func;
}



double init_fermi(int i, int left)
{
  int i_HE, j;	
  double func, den, arg;
	
  if (left)
    j=i+1;
  else
    j=i;
  func=0.;	
	
  for (i_HE=0;i_HE<HIST_SIZES[j];++i_HE){
    arg=(BETAS[i]-BETAS[i+1])*HIST[j][i_HE];
    den=NORM_HIST[j];		
    if (left)
      den*=exp(arg);
    else 
      den*=exp(-arg);
    func+=1./den;		
  }
  if (left)
    return -log(func);
  else
    return log(func);
}
