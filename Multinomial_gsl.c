/*==========================================================
 * Multinomial_gsl.c
 *
 * Return multinomially distributed integers based on input probability vector
 * 
 *========================================================*/
/*  */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h> /* For EXIT_FAILURE, EXIT_SUCCESS */
#include "mat.h"
// #include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#define pi acos(-1.0)




int rlxd_seed(){
    int seed;
    FILE *file_seeds;
    file_seeds=fopen("/dev/urandom","rb"); // /dev/urandom is designed so that each time it is opened it produces a different random seed 
    fread(&seed, sizeof(int), 1, file_seeds);
    fclose(file_seeds);

    return abs(seed);
}


void Multinomial(unsigned int *n, double N, double *x, mwSize K)
{
    
//     Multinomial(n,N,x,K);
    
    
    //    /* calculate eigenvalues and eigenvectors */
    //mexCallMATLAB(2, lhs, 1, &x, "eig");
    
//    printf("N =%f\n",N);   
    
//    printf("K =%d\n",K);      


    
      
    
    //Initialisation for using gls multinomial function
    
    const gsl_rng_type *Type;         //initialise a constant Type which holds type of random number generator
    gsl_rng *r;                              //Initialise pointer to r which is an "instance" of random number generator

    //Set seed for random numbers from current time
//     srand(time(0));
//     srand(rlxd_seed());
    unsigned long int seed = rlxd_seed();         //initialise integer seed and set to value shown
//     unsigned long int seed = random_seed();         //initialise integer seed and set to value shown

    
    gsl_rng_env_setup();                        //read in environmental variables set at command line

    Type = gsl_rng_default;                        //set random number generator to default, which is mt19937
    r = gsl_rng_alloc (Type);                       //Allocate memory for type T generator and get pointer
    gsl_rng_set (r, seed);                      //initialise random number generator

    unsigned int NN = N;                    //Sample size
    size_t KK =K;
//     double mx[KK];                    //Initialise pointer to mx real positive array of probabilities

//     unsigned int n[KK];           //Initialise output integer array
    
//     printf ("generator type: %s\n", gsl_rng_name (r));
//     printf ("seed = %d\n", seed);



//     printf("NN=%d\n",NN); 
//     printf("KK=%d\n",KK); 
    
    
    gsl_ran_multinomial(r, KK, NN, x, n);

      
    gsl_rng_free(r);
    
//     free(n);
//     free(mx);
    
 }   




        
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    
    
     /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:Multinomial_gsl:nrhs","Three inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:Multinomial_gsl:nlhs","One output required.");
    }
    
    /* get input variables */
    double N=mxGetScalar(prhs[0]);
    double *x = mxGetPr(prhs[1]);
    double K=mxGetScalar(prhs[2]);
    
    int *n;
    
    mwSize Ka=K;
    
//     plhs[0] = mxCreateDoubleMatrix(1, Ka,mxREAL);    
    //mxCreateNumericMatrix is used here as we want to create a vector of integers
    plhs[0] = mxCreateNumericMatrix(1, Ka, mxUINT32_CLASS, mxREAL);  


   /* get a pointer to the real data in the output matrix */   
   //As arrays is not DOUBLE then we can't use mxGetPr
   n = (int*) mxGetData(plhs[0]);

   Multinomial(n,N,x,Ka);
    


    
}     