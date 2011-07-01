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
#include "mmap.h"

extern char *get_string(const mxArray *prhs);
extern int mmap_file(const char *filename, int open_mode, void **map, off_t *size) ;
extern int find_interval(unsigned *pos_map, off_t num_entries,	unsigned begin, unsigned end, unsigned *lindex, unsigned *rindex);

/*enum {E_DUMMY, E_SYNTAX};*/
enum {E_SYNTAX = 1, E_BAD_INTERVAL};

/*
 * rhs[0] basename
 *
 * return:
 * lhs[0]: number of stored positions
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	char *basename = get_string(prhs[0]);
	int ret=0 ;
	unsigned *pos_map = NULL;
	off_t pos_size, num_positions=0;
	char filename[MAXLINE], *err_txt;

	if (nlhs!=1 || nrhs<1) {
		ret = -E_SYNTAX;
		err_txt = "fatal: wrong number of arguments\nusage: num_pos = get_num_pos(fname[, interval_list])";
		goto out;
	}


	if (strlen(basename) > MAXLINE - 5){	
		err_txt = "basename too long";
		goto out;
	}

	sprintf(filename, "%s.pos", basename);
	ret = mmap_file(filename, O_RDONLY, (void **)&pos_map, &pos_size);
	err_txt = "pos mmap error";
	if (ret < 0)
		goto out;
	num_positions =	pos_size / sizeof(unsigned);
	/*fprintf(stderr, "num_pos=%i\n", num_positions) ;*/

	if (nrhs>=2 && num_positions>0) {
		if (mxGetM(prhs[1]) != 2) {
			ret = -E_SYNTAX;
			err_txt = "expected N x 2 Matrix\n";
			goto out;
		}

		int total_num_pos = 0;
		double *interval_ptr = mxGetPr(prhs[1]);
		const int num_intervals = mxGetN(prhs[1]);
		int j ; 
		for (j = 0; j < num_intervals; j++) {
			unsigned left = (unsigned)interval_ptr[2 * j];
			unsigned right = (unsigned)interval_ptr[2 * j + 1];
			unsigned lindex, rindex;
			unsigned diff;

			if (left > right) {
				ret = -E_SYNTAX;
				err_txt = "bad interval";
				goto out;
			}
			int num_found = find_interval(pos_map, num_positions, left, right, &lindex, &rindex);
			if (num_found==0)
				continue ;

			/*fprintf(stderr, "left=%i, right=%i, l_idx=%i, r_idx=%i\n", left, right, lindex,rindex) ;*/
			diff = rindex - lindex + 1;
			if (diff<0) {
				fprintf(stderr, "total_num_pos=%i diff=%i  lindex=%i  rindex=%i  num_positions=%i\n", total_num_pos, diff, lindex, rindex, (int) num_positions) ;
				err_txt = "pos mmap error: lindex > rindex!";
				mexErrMsgTxt(err_txt);
			}
			total_num_pos += diff;
		}
		num_positions = total_num_pos ;
	}
	/*fprintf(stderr, "num_pos=%i\n", num_positions) ;*/

	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double * p_num_positions = mxGetPr(plhs[0]) ;
	*p_num_positions = (double) num_positions ;

out:
	if (pos_map)
		munmap(pos_map, pos_size);
	if (ret < 0)
		mexErrMsgTxt(err_txt);
}
