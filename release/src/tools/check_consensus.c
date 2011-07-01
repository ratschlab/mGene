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
/*#define DEBUG 1 */

extern char *get_string(const mxArray *prhs);

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
	char *dna_seq = get_string(prhs[0]);
	int dna_seq_len = mxGetN(prhs[0]);
	
	/* Argument 1*/
	int num_motifs = mxGetNumberOfElements(prhs[1]);
	if (num_motifs==0)
		mexErrMsgTxt("no motives given\n");

	unsigned i=0, j=0, k=0;
	unsigned *motif_length = malloc(num_motifs*sizeof(unsigned));
	if (!motif_length)
		mexErrMsgTxt("ERROR ALLOCATING MEM\n"); 

	char **motifs = NULL;
	motifs = malloc(num_motifs * sizeof(char *));
	if (!motifs)
		mexErrMsgTxt("ERROR ALLOCATING MEM\n");	

	for (i = 0; i < num_motifs; i++)
	{
		motifs[i] = get_string(mxGetCell(prhs[1], i));
		j = 0;
		while (motifs[i][j] != 0)
			j++;
		motif_length[i] = j;
	}
	
	/* Argument 2*/
	int num_offsets = mxGetN(prhs[2]);
	assert(num_offsets==num_motifs);
	double *offsets = mxGetPr(prhs[2]);  

	/* Argument 3*/
	int num_cands = mxGetN(prhs[3]);
	int M_cands = mxGetM(prhs[3]);
	double *cands = mxGetPr(prhs[3]);  

	if (num_cands==0||M_cands==0)/*0x1 or 1x0 vector*/
		num_cands = 0;
#ifdef DEBUG
	printf("nlhs: %i, nrhs: %i\n",nlhs, nrhs);
	printf("num_cands:%i\n",num_cands);
	printf("dna_seq_len: %i\n",dna_seq_len);
	printf("num_motivs: %i\n",num_motifs);
	printf("num_offsets:%i\n",num_offsets);
	for (i = 0; i < num_motifs; i++)
	{
		printf("motif %i: %s\n", i, motifs[i]);
		printf("motif_length[%i]:%i\n",i, motif_length[i]); 
	}
#endif

	/* find true positions */
	double *true_pos = malloc(num_cands*sizeof(double));
	if (!true_pos)
		mexErrMsgTxt("ERROR ALLOCATING MEM\n");
	for (i=0; i<num_cands; i++)
	{
		true_pos[i] = 0;
		for (j=0; j<num_motifs; j++)
		{
			if (true_pos[i] == 1)
				break;
			for (k=0; k<motif_length[j]; k++)
			{
				if (cands[i]+k-offsets[j]<0)
					break;
				int pos = (int) cands[i]+k- (int) offsets[j];
				if (pos>=dna_seq_len)
				{
					/*printf("check_consensus.c: position (%i) larger than length of the sequence (%i); this position will be removed\n", pos, dna_seq_len);*/
					true_pos[i] = 0;
					break;
				}
				if (dna_seq[pos]==motifs[j][k])
					true_pos[i] = 1;
				else
				{
					true_pos[i] = 0;
					break;	
				}
			}
		}
	}
	/* handle output*/
	plhs[0] = mxCreateDoubleMatrix(1, num_cands, mxREAL);
	double *out_ptr = mxGetPr(plhs[0]);

	for (i=0; i<num_cands; i++)
		out_ptr[i] = true_pos[i];

	free(motif_length);
	for (j=0; j<num_motifs; j++)
		free(motifs[j]);
	free(motifs);
	free(true_pos);

}

