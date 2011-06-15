/*
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Written (W) 2005-2009 Gunnar Raetsch, Gabriele Schweikert, Jonas Behr, Soeren Sonnenburg, Alexander Zien, Georg Zeller, Andre Noll, Cheng Soon Ong, Petra Philips
  Copyright (C) 2005-2009 Max Planck Society and Fraunhofer Institute FIRST 
*/

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

extern char *get_string(const mxArray *prhs);
extern int mmap_file(const char *filename, int open_mode, void **map, off_t *size) ;
extern int find_interval(unsigned *pos_map, off_t num_entries,	unsigned begin, unsigned end, unsigned *lindex, unsigned *rindex);

/*enum {E_DUMMY, E_SYNTAX};*/
enum {E_SYNTAX = 1, E_BAD_INTERVAL};

/*
 * rhs[0] basename
 * rhs[1] list of m score names
 * rsh[2] Interval Matrix (n x 2)
 *
 * return:
 * lhs[0]: list of positions that lie in any given interval, k := number of such positions
 * lhs[1]: k x m score Matrix
 * lhs[2]: index vector (k x 1), lhs[2][j] = l <=> rhs[2][l, 0] <= lhs[0][j] <= rhs[2][l, 1].
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	char filename[MAXLINE], *err_txt;
	int ret;
	char *basename = get_string(prhs[0]);
	unsigned i, j, num_scores = 0;
	float **score_maps = NULL;
	off_t *score_sizes = NULL;
	char **score_names = NULL;
	double *interval_ptr;
	const int num_intervals = mxGetN(prhs[2]);
	unsigned m, k,total_num_pos, *pos_map = NULL;
	off_t pos_size, num_positions;
	double *pos_ptr, *score_ptr, *index_ptr;

	if ((nlhs<2) || (nrhs!=3))
	  {
	    ret = -E_SYNTAX;
	    err_txt = "fatal: wrong number of arguments\nusage: [pos, scores, idx] = get_num_pos(fname, score_names, positions)";
	    goto out;
	  }

	if (!mxIsCell(prhs[1]))
	  {
	    ret = -E_SYNTAX;
	    err_txt = "fatal:expected a cell as second argument (score_names)\nusage: [pos, scores, idx] = get_num_pos(fname, score_names, positions)";
	    goto out;
	  }

	ret = -E_SYNTAX;
	err_txt = "fatal: no interval(s) given";
	if (!num_intervals)
		goto out;
	err_txt = "fatal: no score file(s) given";
	num_scores = mxGetNumberOfElements(prhs[1]);
	if (!num_scores)
		goto out;
	score_maps = malloc(num_scores * sizeof(void *));
	memset(score_maps, 0, num_scores * sizeof(void *));
	score_sizes = malloc(num_scores * sizeof(off_t *));
	score_names = malloc(num_scores * sizeof(char *));
	memset(score_names, 0, num_scores * sizeof(char *));
	for (i = 0; i < num_scores; i++)
		score_names[i] = get_string(mxGetCell(prhs[1], i));
	err_txt = "expected N x 2 Matrix\n";
	if (mxGetM(prhs[2]) != 2)
		goto out;
	sprintf(filename, "%s.pos", basename);
	ret = mmap_file(filename, O_RDONLY, (void **)&pos_map, &pos_size);
	err_txt = "pos mmap error";
	if (ret < 0)
		goto out;
	num_positions =	pos_size / sizeof(unsigned);
	
#ifdef DEBUG
	err_txt = "pos map is not sorted";
	printf("num_positions: %i\n", num_positions);
	for (i=1; i<num_positions; ++i) {
	    if (pos_map[i-1] > pos_map[i])
	            goto out;
	}
	printf("successfully checked that pos_map is sorted\n");
#endif
	
	for (i = 0; i < num_scores; i++) {
		sprintf(filename, "%s.%s", basename, score_names[i]);
		ret = mmap_file(filename, O_RDONLY, (void **)&score_maps[i], &score_sizes[i]);
		err_txt = "score mmap error";
		if (ret < 0)
		        goto out;
	}
	total_num_pos = 0;
	interval_ptr = mxGetPr(prhs[2]);
	for (j = 0; j < num_intervals; j++) {
	    unsigned left = (unsigned)interval_ptr[2 * j];
	    unsigned right = (unsigned)interval_ptr[2 * j + 1];
#ifdef DEBUG
	      for (i=0; i<num_positions; ++i) {
		if (left <= pos_map[i] && pos_map[i] <= right)
		  printf("%i ", pos_map[i]);
	      }
	      printf("\n\n");
#endif
	    unsigned lindex, rindex;
	    
	    ret = -E_SYNTAX;
	    err_txt = "bad interval";
	    if (left > right)
	        goto out;
	    ret = find_interval(pos_map, num_positions, left, right, &lindex, &rindex);
	    if (ret > 0) {/* non-empty interval */
	        total_num_pos += rindex - lindex + 1;
#ifdef DEBUG
		printf("lidx: %u, ridx: %u, t: %u\n", lindex, rindex, total_num_pos);
#endif
	    }
	}
	plhs[0] = mxCreateDoubleMatrix(total_num_pos, 1, mxREAL);
	pos_ptr = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(total_num_pos, num_scores, mxREAL);
	score_ptr = mxGetPr(plhs[1]);
	if (nlhs>=2) {
	    plhs[2] = mxCreateDoubleMatrix(total_num_pos, 1, mxREAL);
	    index_ptr = mxGetPr(plhs[2]);
	} else
	    index_ptr = NULL;
	k = 0;
	for (j = 0; j < num_intervals; j++) {
	    unsigned left = (unsigned)interval_ptr[2 * j];
	    unsigned right = (unsigned)interval_ptr[2 * j + 1];
	    unsigned lindex, rindex;
	    
	    ret = find_interval(pos_map, num_positions, left, right, &lindex, &rindex);
	    if (ret <= 0) /* empty interval */
		continue;
            assert(ret == rindex - lindex + 1);
	    for (m = 0; m < ret; m++) {
	      if (index_ptr)
	  	index_ptr[k + m] = j;
	      err_txt = "returned position out of interval";
              if (left > pos_map[lindex + m] || pos_map[lindex + m] > right)
		goto out;

#ifdef DEBUG
              printf("l: %u, r: %u, p: %u\n", left, right, pos_map[lindex + m]);
#endif
	      pos_ptr[k + m] = pos_map[lindex + m];
	      for (i = 0; i < num_scores; i++)
		/*score_ptr[i + num_scores * (k + m)] = (double)score_maps[i][lindex + m];*/
		score_ptr[i * total_num_pos + (k + m)] = (double)score_maps[i][lindex + m];
	    }
	    k += ret;
	}
	ret = 1;
 out:
	for (i = 0; i < num_scores; i++) {
	  if (score_maps[i])
	    munmap(score_maps[i], score_sizes[i]);
	  free(score_names[i]);
	}
	free(score_maps);
	free(score_names);
	free(score_sizes);
	if (pos_map)
	  munmap(pos_map, pos_size);
	if (ret < 0)
	  mexErrMsgTxt(err_txt);
}
