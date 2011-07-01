MEX_BIN ?= /fml/ag-raetsch/share/software/matlab-7.6/bin/mex
MATLAB_INCDIR ?= /fml/ag-raetsch/share/software/matlab/extern/include
MKOCT_BIN ?= mkoctfile
CPLEX_INCL ?= /fml/ag-raetsch/share/software/ilog/cplex91/include/ilcplex/
CPLEX_LIB ?= /fml/ag-raetsch/share/software/ilog/cplex91/lib/x86-64_sles9.0_3.3/static_pic/

all: qp_solve.mexa64 qp_solve.mex cplex_init_quit.mexa64 cplex_init_quit.mex

qp_solve.mexa64:
	${MEX_BIN} -g -O qp_solve.c -I$(MATLAB_INCDIR) -I$(CPLEX_INCL) -L$(CPLEX_LIB) -lcplex

cplex_init_quit.mexa64: cplex_init_quit.c
	${MEX_BIN} -g -O cplex_init_quit.c -I$(MATLAB_INCDIR) -I$(CPLEX_INCL) -L$(CPLEX_LIB) -lcplex

qp_solve.mex:
	${MKOCT_BIN} -g --mex qp_solve.c -I$(CPLEX_INCL) -L$(CPLEX_LIB) -lcplex

cplex_init_quit.mex: cplex_init_quit.c
	${MKOCT_BIN} -g --mex cplex_init_quit.c -I$(CPLEX_INCL) -L$(CPLEX_LIB) -lcplex

