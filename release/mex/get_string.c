/*
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Written (W) 2005-2009 Gunnar Raetsch, Gabriele Schweikert, Jonas Behr, Soeren Sonnenburg, Alexander Zien, Georg Zeller, Andre Noll, Cheng Soon Ong, Petra Philips
  Copyright (C) 2005-2009 Max Planck Society and Fraunhofer Institute FIRST 
*/

#include <mex.h>
#include <stdlib.h>
char *get_string(const mxArray *prhs)
{
	char *buf;
	int buflen;

	if (!prhs)
		mexErrMsgTxt("get_string called with NULL pointer arg");
	if (!mxIsChar(prhs))
		mexErrMsgTxt("input is not a string");
	if (mxGetM(prhs) != 1)
		mexErrMsgTxt("input is not a row vector");
	buflen = mxGetN(prhs) + 1;
	buf = malloc(buflen);
	/* copy the string from prhs into buf and add terminating NULL char */
	if (mxGetString(prhs, buf, buflen))
		mexErrMsgTxt("not enough space");
	return buf;
}
