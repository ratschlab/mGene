all: make_subfolders check_solver check_shogun

include mgene_config.sh
export GENOME="../../genome"
export SAMDIR="../../samtools"

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
ifeq ($(OPTIMIZER),cplex)
	@echo "##################################################" 
	@echo "# compiling QP solver tools"
	@echo "##################################################" 
	$(MAKE) -C solver/
endif

	
## check shogun
################################################################################
check_shogun:
	@echo "##################################################" 
	@echo "# Trying to run shogun binary ..."
	@echo "# $(shell sh check_shogun.sh)!"
	@echo  "##################################################" 

## check solver
################################################################################
check_solver:
	@echo "##################################################" 
	@echo "# Trying to run QP solver ..."
	@echo "# $(shell sh check_solver.sh) !"
	@echo  "##################################################" 

