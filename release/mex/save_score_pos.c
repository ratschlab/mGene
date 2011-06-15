/*
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Written (W) 2005-2009 Gunnar Raetsch, Gabriele Schweikert, Jonas Behr, Soeren Sonnenburg, Alexander Zien, Georg Zeller, Andre Noll, Cheng Soon Ong, Petra Philips
  Copyright (C) 2005-2009 Max Planck Society and Fraunhofer Institute FIRST 
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <errno.h>

#include <mex.h>

#define MAXLINE 255

/* define this empty to get rid of log messages */
#define LOG printf

extern char *get_string(const mxArray *prhs);

enum {E_SYNTAX = 1, E_OPEN, E_WRITE, E_BAD_FILENAME, E_MONOTONICITY};
/*
 * We fundamentally don't like some basenames: We don't like
 * dot or dot-dot, and for obvious reasons don't want any
 * slashes.
 */
int verify_basename(const char *name)
{
	if (!name)
		return -E_BAD_FILENAME;
	if (!*name)
		return -E_BAD_FILENAME;
	if (!strcmp(name, "."))
		return -E_BAD_FILENAME;
	if (!strcmp(name, ".."))
		return -E_BAD_FILENAME;
	if (strchr(name, '/'))
		return -E_BAD_FILENAME;
	return 1;
}

/**
 * compute the difference of two time values
 *
 * \param b minuend
 * \param a subtrahend
 * \param diff difference is stored here
 *
 * \return If b > a, compute b - a and return 1. if b < a
 * compute a - b and return -1. If diff is not  NULL, diff gets
 * filled in with the absolute value |b - a|.
 */
int tv_diff(const struct timeval *b, const struct timeval *a, struct timeval *diff)
{
	int ret = 1;

	if ((b->tv_sec < a->tv_sec) ||
			((b->tv_sec == a->tv_sec) && (b->tv_usec < a->tv_usec))) {
		const struct timeval *tmp = a;
		a = b;
		b = tmp;
		ret = -1;
	}
	if (!diff)
		return ret;
	diff->tv_sec = b->tv_sec - a->tv_sec;
	if (b->tv_usec < a->tv_usec) {
		diff->tv_sec--;
		diff->tv_usec = 1000 * 1000 - a->tv_usec + b->tv_usec;
	} else
		diff->tv_usec = b->tv_usec - a->tv_usec;
	return ret;
}

long unsigned tv2ms(const struct timeval *tv)
{
	return tv->tv_sec * 1000 + (tv->tv_usec + 500) / 1000;
}



/*
 * n: number of positions
 * m: number of value columns (scores)
 * rhs[0]: vector of n positions
 * rhs[1]: n x m matrix of score values
 * rhs[2]: basename
 * rhs[3]: vector containing m filenames for the output files
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n;
	char filename[MAXLINE], *path_prefix = NULL, **score_names = NULL, *err_txt = NULL;
	int i, ret, pos_fd = -1, m = -1;
	int *score_fds = NULL;
	double *pos_ptr, *score_ptr;
	struct timeval start, end, diff;

	/*LOG("%s startup\n", __func__);*/
	gettimeofday(&start, NULL);
	ret = -E_SYNTAX;
	err_txt = "expected 4 arguments\n";
	if (nrhs != 4)
		goto out;
	n = mxGetN(prhs[0]);
	LOG("num_positions: %d, num_scores: %d\n", n, mxGetN(prhs[1]));
	err_txt = "fatal: num_positions != num_scores\n";
	if (n != mxGetN(prhs[1]))
		goto out;
	m = mxGetM(prhs[1]);
	LOG("num score value colums: %d, num oputput files: %d\n", m, mxGetN(prhs[3]));
	err_txt = "fatal: num score value columns != num output files\n";
	if (mxGetN(prhs[3]) != m && n!=0)
		goto out;
    /* if n=0, then the m-dimensionality is often not informative ... */
	m=mxGetN(prhs[3]) ;
	
	path_prefix = get_string(prhs[2]);
	LOG("path_prefix: %s\n", path_prefix);
	score_names = malloc(m * sizeof(char*));
	memset(score_names, 0, m * sizeof(char*));
	score_fds = malloc(m * sizeof(int));
	LOG("getting %d score names\n", m);
	for (i = 0; i < m; i++) {
		score_names[i] = get_string(mxGetCell(prhs[3], i));
		LOG("score name %d: %s\n", i, score_names[i]);
		ret = verify_basename(score_names[i]);
		err_txt = "invalid score name";
		if (ret < 0)
			goto out;
		score_fds[i] = -1;
	}
	sprintf(filename, "%s.pos", path_prefix);
	ret = open(filename, O_WRONLY | O_CREAT | O_EXCL, 0666);
	if (ret < 0) {
	  if (errno == EEXIST)
		err_txt = "position output file already exists\n";
	  else
		err_txt = "could not open position output file\n";
	  ret = -E_OPEN;
	  goto out;
	}
	pos_fd = ret;
	for (i = 0; i < m; i++) {
		sprintf(filename, "%s.%s", path_prefix, score_names[i]);
		ret = open(filename, O_WRONLY | O_CREAT | O_EXCL, 0666);
		if (ret < 0) {
		  if (errno == EEXIST)
		    err_txt = "score output file already exists\n";
		  else
		    err_txt = "could not open score output file\n";
		  ret = -E_OPEN;
		  goto out;
		}
		score_fds[i] = ret;
	}
	pos_ptr = mxGetPr(prhs[0]);
	score_ptr = mxGetPr(prhs[1]);
	LOG("writing pos file and %d score files, each %d bytes...\n",
		m, n * sizeof(float));

        unsigned *pos_buf = (unsigned*) malloc(n*sizeof(unsigned)) ;
	for (i = 0; i < n; i++) {
		unsigned pos = (unsigned) pos_ptr[i];
		if (i && pos_ptr[i] < pos_ptr[i - 1]) {
			ret = -E_MONOTONICITY;
			LOG("pos[%d] = %u < %u = pos[%d]", i, pos_ptr[i],
				pos_ptr[i - 1], i - 1);
			err_txt = "positions are increasing!\n";
			goto out;
		}
                pos_buf[i] = pos ;
	}
	ret = write(pos_fd, pos_buf, sizeof(unsigned)*n);
	free(pos_buf) ;
	if (ret < 0) {
		ret = -E_WRITE;
		err_txt = "pos file write error\n";
		goto out;
	}
        float *score_buf = (float*) malloc(n*sizeof(float)) ;
	int j ;
	for (j = 0; j < m; j++) {
	  for (i = 0; i < n; i++) {
	    float score = (float) score_ptr[i * m + j];
	    score_buf[i]=score ;
	  }
	  ret = write(score_fds[j], score_buf, sizeof(float)*n);
	  if (ret < 0) {
	    ret = -E_WRITE;
	    err_txt = "score file write error\n";
	    free(score_buf) ;
	    goto out;
	  }
	}
	free(score_buf) ;

	ret = 1;
out:
	if (score_names) {
		for (i = 0; i < m; i++) {
			free(score_names[i]);
			if (score_fds[i] >= 0)
				close(score_fds[i]);
		}
		free(score_names);
	}
	free(score_fds);
	free(path_prefix);
	if (pos_fd >= 0)
		close(pos_fd);
	gettimeofday(&end, NULL);
	tv_diff(&end, &start, &diff);
	if (ret < 0)
		mexErrMsgTxt(err_txt);
	LOG("success (took %lums)\n", tv2ms(&diff));
}
