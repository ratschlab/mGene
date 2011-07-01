#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/mman.h>
#include <time.h>
#include <mex.h>
#include <assert.h>

#define MAXLINE 255
/*#define DEBUG 1*/ 


/*
 * prhs[0] DNA sequence
 * prhs[1] cell of motivs
 * prhs[2] offsets of motivs
 * prhs[3] candidate positions
 *
 * return:
 * plhs[0] validation vector
 *
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	/* Argument 0*/
	int num_values = mxGetN(prhs[0]);
	double *values = mxGetPr(prhs[0]);  

	/* Argument 1*/
	int num_prob = mxGetN(prhs[1]);
	double *prob = mxGetPr(prhs[1]);  

	/* Argument 2*/
	int num_limits = mxGetN(prhs[2]);
	double *limits = mxGetPr(prhs[2]);  



#ifdef DEBUG
	printf("nlhs: %i, nrhs: %i\n",nlhs, nrhs);
	printf("num_values:%i\n",num_values);
	printf("num_prob: %i\n",num_prob);
	printf("num_limits: %i\n",num_limits);
	/*printf("num_offsets:%i\n",num_offsets);
	for (i = 0; i < num_motifs; i++)
	{
		printf("motif %i: %s\n", i, motifs[i]);
		printf("motif_length[%i]:%i\n",i, motif_length[i]); 
	}*/
#endif
	/*check input data 
	 *************************/
	if (num_values==0)
		mexErrMsgTxt("Error: array of values for transformation is empty\n");
	if (num_limits<1)
		mexErrMsgTxt("Error: no support points for transformation given\n");
	if (num_prob<1)
		mexErrMsgTxt("Error: no plif values for transformation given\n");
        if (num_limits!=num_prob)
	{
		fprintf(stdout, "num_limits:%i num_prob:%i \n", num_limits, num_prob);
		mexErrMsgTxt("Error: number of plif values does not match number of support points\n");
	}
	double val;
        int idx;
	plhs[0] = mxCreateDoubleMatrix(1, num_values, mxREAL);
        double* output = mxGetPr(plhs[0]);
	if (output == NULL)
	{
		mexErrMsgTxt("OUT OF MEMORY");
		exit(EXIT_FAILURE);
	}
	int i;
    	for(i=0; i<num_values;i++)
	{
		val = values[i];
		idx = 0;
		int j;
		for (j=0; j<num_limits; j++)
		{
        	        if (limits[j]<=val)
                	        idx++ ;
                	else
                        	break ; /* assume it is monotonically increasing*/
         	}
	        if (idx==0)
	                output[i]=prob[0] ;
	        else if (idx==num_limits)
	                output[i]=prob[num_limits-1] ;
	        else
	        {
	                output[i] = (prob[idx]*(val-limits[idx-1]) + prob[idx-1]*
	                           (limits[idx]-val)) / (limits[idx]-limits[idx-1]) ;  
	        }
	}
}

