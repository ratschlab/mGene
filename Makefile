all: make_subfolders check_solver check_shogun

include mgene_config.sh
export GENOME="../../genome"
export SAMDIR="../../samtools"

RES_SG := $(shell sh check_shogun.sh)
RES_QP := $(shell sh check_solver.sh)


## compile subfolders
################################################################################
make_subfolders:
	@echo "##################################################" 
	@echo "# compiling ngs tools"
	@echo "##################################################" 
	$(MAKE) -C shared_tools/ngs/
	@echo "##################################################" 
	@echo "# mGene tools and utils"
	@echo "##################################################" 
	$(MAKE) -C tools/
	@echo "##################################################" 
	@echo "# compiling ngs tools"
	@echo "##################################################" 
	$(MAKE) -C samtools/
	@echo "##################################################" 
	@echo "# compiling QP solver tools"
	@echo "##################################################" 
	$(MAKE) -C solver/

	
## check shogun
################################################################################
check_shogun:
	@echo "##################################################" 
	@echo "# Trying to run shogun binary ..."
ifeq ($(RES_SG),sg_test_succ)
	@echo "# success!"
else
	@echo "# failure!"
endif
	@echo  "##################################################" 

## check solver
################################################################################
check_solver:
	@echo "##################################################" 
	@echo "# Trying to run QP solver ..."
	$(if $(eq $(shell sh check_solver.sh),"qp_test_succ"),@echo "# success!",@echo $(shell sh check_solver.sh))
ifeq ($(shell sh check_solver.sh),qp_test_succ)
	@echo "# success!"
else
	@echo "# failure!"
endif
	@echo  "##################################################" 

