
BAM_ORIG=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans/polyA_left_sam_stranded.mapped.2.bam
BAM_SMALL=~/svn/releases/mGeneToolbox-0.2.0/release/examples/data/elegans_chrII.sam
BAM_DIR=~/svn/releases/mGeneToolbox-0.2.0/release/examples/data/
FN_FASTA=/fml/ag-raetsch/share/projects/rgasp/genomes/elegans/c_elegans.WS200.dna.fa

samtools view $BAM_ORIG II:1-150000 > $BAM_SMALL

~/svn/tools/ngs/sam_to_bam.sh $BAM_DIR $FN_FASTA

rm release/examples/data/elegans_chrII.sam.gz
