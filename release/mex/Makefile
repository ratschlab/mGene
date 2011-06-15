mex=mex
mkoctfile=mkoctfile 

all: save_score_pos.mexa64 interval_query.mexa64 get_num_pos.mexa64 check_consensus.mexa64 save_score_pos.mex interval_query.mex get_num_pos.mex check_consensus.mex

octave: save_score_pos.mex interval_query.mex get_num_pos.mex check_consensus.mex
matlab: save_score_pos.mexa64 interval_query.mexa64 get_num_pos.mexa64 check_consensus.mexa64 

save_score_pos.mexa64: save_score_pos.c get_string.c 
	rm -f *.o
	${mex} -g -O save_score_pos.c get_string.c -I/fml/ag-raetsch/share/software/matlab/extern/include
interval_query.mexa64: get_string.c interval_query.c mmap.c find_interval.c
	rm -f *.o
	${mex} -g -O interval_query.c get_string.c find_interval.c mmap.c -I/fml/ag-raetsch/share/software/matlab/extern/include
get_num_pos.mexa64: get_string.c get_num_pos.c mmap.c find_interval.c
	rm -f *.o
	${mex} -g -O get_num_pos.c get_string.c find_interval.c mmap.c -I/fml/ag-raetsch/share/software/matlab/extern/include
check_consensus.mexa64: get_string.c check_consensus.c
	rm -f *.o
	${mex} -g -O check_consensus.c get_string.c -I/fml/ag-raetsch/share/software/matlab/extern/include

save_score_pos.mex: save_score_pos.c get_string.c 
	rm -f *.o
	mkoctfile -g --mex save_score_pos.c get_string.c 
interval_query.mex: get_string.c interval_query.c mmap.c find_interval.c
	rm -f *.o
	mkoctfile -g --mex interval_query.c get_string.c find_interval.c mmap.c 
get_num_pos.mex: get_string.c get_num_pos.c mmap.c find_interval.c
	rm -f *.o
	mkoctfile -g --mex get_num_pos.c get_string.c find_interval.c mmap.c 
check_consensus.mex: get_string.c check_consensus.c
	rm -f *.o
	mkoctfile -g --mex check_consensus.c get_string.c

