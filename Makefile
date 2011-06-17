all: check_shogun make_subfolders

include mgene_config.sh
export GENOME="../../genome"
export SAMDIR="../../samtools"

RES := $(shell sh check_shogun.sh)


## compile subfolders
################################################################################
make_subfolders:
	@echo
	@echo compiling ngs tools
	@echo
	$(MAKE) -C shared_tools/ngs/
	@echo
	@echo mGene tools and utils
	@echo
	$(MAKE) -C tools/
	@echo
	@echo compiling ngs tools
	@echo
	$(MAKE) -C samtools/
	
## check shogun
################################################################################
check_shogun:
	@echo Checking shogun binary: ...
ifeq ($(RES),sg_test_succ)
	@echo success
else
	@echo failed
endif

## check solver
################################################################################
