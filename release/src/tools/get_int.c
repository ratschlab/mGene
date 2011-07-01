#include <mex.h>
#include <stdlib.h>

#include <stdio.h>

int get_int(const mxArray *prhs)
{
	if (!prhs)
		mexErrMsgTxt("get_int: called with NULL pointer arg");
	if (!mxIsDouble(prhs))
		mexErrMsgTxt("get_int: input is not a string");
	if (mxGetN(prhs)!=1 || mxGetM(prhs)!=1)
		 mexErrMsgTxt( "get_int: argument should be a scalar\n");

	double* arg = mxGetPr(prhs);	
	return (int) *arg;
}
