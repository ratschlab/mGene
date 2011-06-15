
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
rsync -a ~/svn/projects/splicing/splicegraphs release/src/
rsync -a ~/svn/tools/ngs release/src/shared_tools/
rsync -a ~/svn/tools/genomes release/src/shared_tools/
rsync -a ~/svn/tools/rproc release/src/shared_tools/
rsync -a ~/svn/tools/utils release/src/shared_tools/
rsync -a ~/svn/tools/prob_plif release/src/shared_tools/

cp ~/svn/projects/mGene_core/auxiliary_data/add_reads_from_bam.m release/src/auxiliary_data/

# modified files
cp get_base_dir.m_template release/src/get_base_dir.m
cp paths.m_template release/src/paths.m
cp configure release/src/configure
