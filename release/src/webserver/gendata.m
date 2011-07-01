function gendata(organism, fasta_fname, gff3_fname, info_fname)
% galaxy_gendata(organism, fasta_fname, gff3_fname, info_fname)

global MGENE_GALAXY_DIR

unix(sprintf('cp %s/examples/%s.fasta %s', MGENE_GALAXY_DIR, organism, fasta_fname)) ;
unix(sprintf('cp %s/examples/%s.gff3 %s', MGENE_GALAXY_DIR, organism, gff3_fname)) ;
unix(sprintf('cp %s/examples/%s.info %s', MGENE_GALAXY_DIR, organism, info_fname)) ;
