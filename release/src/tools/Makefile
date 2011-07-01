MEX_BIN ?= /fml/ag-raetsch/share/software/matlab-7.6/bin/mex
MATLAB_INCDIR ?= /fml/ag-raetsch/share/software/matlab/extern/include
MKOCT_BIN ?= mkoctfile
mex_incl = "-I$(MATLAB_INCDIR)"
all: unique_transcripts.mexa64 interval_overlap.mexa64 plif_transform.mexa64 save_score_pos.mexa64 interval_query.mexa64 get_num_pos.mexa64 check_consensus.mexa64 unique_transcripts.mex interval_overlap.mex plif_transform.mex save_score_pos.mex interval_query.mex get_num_pos.mex check_consensus.mex 

unique_transcripts.mexa64: unique_transcripts.cpp
	${MEX_BIN} -g -O unique_transcripts.cpp ${mex_incl}
interval_overlap.mexa64: interval_overlap.cpp
	${MEX_BIN} -g -O interval_overlap.cpp ${mex_incl}
plif_transform.mexa64: plif_transform.c
	${MEX_BIN} -g -O plif_transform.c ${mex_incl}
save_score_pos.mexa64: save_score_pos.c get_string.c 
	${MEX_BIN} -g -O save_score_pos.c get_string.c ${mex_incl}
interval_query.mexa64: get_string.c interval_query.c mmap.c find_interval.c
	${MEX_BIN} -g -O interval_query.c get_string.c find_interval.c mmap.c ${mex_incl}
get_num_pos.mexa64: get_string.c get_num_pos.c mmap.c find_interval.c
	${MEX_BIN} -g -O get_num_pos.c get_string.c find_interval.c mmap.c ${mex_incl}
check_consensus.mexa64: get_string.c check_consensus.c
	${MEX_BIN} -g -O check_consensus.c get_string.c ${mex_incl}
#
unique_transcripts.mex: unique_transcripts.cpp
	${MKOCT_BIN} -g --mex unique_transcripts.cpp
interval_overlap.mex: interval_overlap.cpp
	${MKOCT_BIN} -g --mex interval_overlap.cpp
plif_transform.mex: plif_transform.c
	${MKOCT_BIN} -g --mex plif_transform.c
save_score_pos.mex: save_score_pos.c get_string.c 
	${MKOCT_BIN} -g --mex save_score_pos.c get_string.c 
interval_query.mex: get_string.c interval_query.c mmap.c find_interval.c
	${MKOCT_BIN} -g --mex interval_query.c get_string.c find_interval.c mmap.c 
get_num_pos.mex: get_string.c get_num_pos.c mmap.c find_interval.c
	${MKOCT_BIN} -g --mex get_num_pos.c get_string.c find_interval.c mmap.c 
check_consensus.mex: get_string.c check_consensus.c
	${MKOCT_BIN} -g --mex check_consensus.c get_string.c

