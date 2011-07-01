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
/* #define DEBUG 1 */

bool isequal(double* trans1, int len1, double* trans2, int len2);

/*
 * prhs[0] cell array of transcripts
 * prhs[1] number of transcripts
 *
 * return:
 * plhs[0] index of unique transcripts
 *
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs!=1)
		mexErrMsgTxt("Expected one arg: cell array of transcripts \n");

	/* Argument 0*/
	int num_transcripts = mxGetNumberOfElements(prhs[0]);
	//printf("num_transcripts: %i\n", num_transcripts);

	if (num_transcripts==0)
	{
		plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
 		return;
	}
	int num_pos[num_transcripts];
	double* transcripts[num_transcripts];

	for (int i = 0; i<num_transcripts; i++)
	{
		num_pos[i] = (int)  mxGetM(mxGetCell(prhs[0], i));
		int len_row = (int)  mxGetN(mxGetCell(prhs[0], i));
		transcripts[i] =  mxGetPr(mxGetCell(prhs[0], i));
		if (len_row<2 || len_row>3)
			mexErrMsgTxt("expected two or three values per row\n");

		
		num_pos[i] *= len_row;
		//printf("%i %i\n",i , num_pos[i]);
	}		
	
	plhs[0] = mxCreateDoubleMatrix(1, num_transcripts, mxREAL);
	double* ret = mxGetPr(plhs[0]);

	ret[0] = 1;
	for (int i=num_transcripts-1;i>0;i--)
	{
		ret[i] = 1;
		for (int j=i-1;j>=0;j--)
			if (isequal(transcripts[i], num_pos[i],transcripts[j], num_pos[j]))
			{
				ret[i] = 0;
			}
	}
}

bool isequal(double* trans1, int len1, double* trans2, int len2)
{
	if (len1!=len2)
		return false;

	for (int i=0; i<len1; i++)
		if (((int)trans1[i])!=((int)trans2[i]))
			return false;
	return true;
}
