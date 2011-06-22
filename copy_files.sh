
rsync -a ~/svn/projects/mGene_core/configuration release/src/
rsync -a ~/svn/projects/mGene_core/contents release/src/
rsync -a ~/svn/projects/mGene_core/data_processing_lsl release/src/
rsync -a ~/svn/projects/mGene_core/data_processing_signals release/src/
rsync -a ~/svn/projects/mGene_core/eval_lsl release/src/
rsync -a ~/svn/projects/mGene_core/genome release/src/
rsync -a ~/svn/projects/mGene_core/lsl release/src/
rsync -a ~/svn/projects/mGene_core/signals release/src/
rsync -a ~/svn/projects/mGene_core/tools release/src/
rsync -a ~/svn/projects/mGene_core/utils release/src/
rsync -a ~/svn/projects/mGene_core/webserver release/src/
rsync -a ~/svn/projects/mGene_core/solver release/src/
rsync -a ~/svn/projects/splicing/splicegraphs release/src/
rsync -a ~/svn/tools/ngs release/src/shared_tools/
rsync -a ~/svn/tools/genomes release/src/shared_tools/
rsync -a ~/svn/tools/rproc release/src/shared_tools/
rsync -a ~/svn/tools/utils release/src/shared_tools/
rsync -a ~/svn/tools/prob_plif release/src/shared_tools/
rsync -a /fml/ag-raetsch/share/software/samtools/ release/src/samtools/

cp ~/svn/projects/mGene_core/auxiliary_data/add_reads_from_bam.m release/src/auxiliary_data/
GFF_DIR2=~/svn/projects/genefinding/parsegff
GFF_DIR=~vipin/svn/projects/genefinding/parsegff
cp $GFF_DIR/CurateGenes.m release/src/parsegff/
cp $GFF_DIR/GFFParser.py release/src/parsegff/
cp $GFF_DIR/arrange_genes.m release/src/parsegff/
cp GFFParser.sh release/src/parsegff/

# modified files
cp get_base_dir.m_template release/src/get_base_dir.m
cp paths.m_template release/src/paths.m
cp shogun_settings.m_template release/src/shogun_settings.m
cp Makefile release/src/Makefile
cp configure release/src/configure
cp cplex_license.m release/src/solver/cplex_license.m
cp check_shogun.sh release/src/check_shogun.sh
cp check_solver.sh release/src/check_solver.sh
cp rproc_policy.m release/src/configuration/rproc_policy.m
