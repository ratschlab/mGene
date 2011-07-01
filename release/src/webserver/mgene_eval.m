function mgene_eval(label_fn, label_dir, version, seqid, source, level, ignore_missing_referenced_entries, run_locally, converter, ...
			predictions_fn, predictions_dir, genome_fasta, performance_fn, performance_dir,  eval_genomewide)

label_fn, label_dir, version, seqid, source, level, ignore_missing_referenced_entries, run_locally, converter,
predictions_fn, predictions_dir, genome_fasta, performance_fn, performance_dir,  eval_genomewide

[ret, timedate] = unix('date') ;
fprintf('started mgene_eval at %s', timedate) ;

fprintf('Results are written to: %s\n', eval_genomewide) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process fasta file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

information_file = [performance_dir '/genome/info.txt'] ;
genome_dir = [performance_dir '/genome/'] ;
fasta_dir = [performance_dir '/fasta/'] ;
nice_mkdir(genome_dir) ;
nice_mkdir(fasta_dir) ;

fprintf('Creating genome information object\n') ;
galaxy_genometool(genome_fasta, fasta_dir, information_file, genome_dir) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% processing gff-file no 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

label_genes_dir = [label_dir '/anno/'];
fn_info = [label_dir '/anno/out'];

galaxy_gff2anno(genome_dir, label_fn, fn_info, label_genes_dir, version, seqid, source, level, ignore_missing_referenced_entries, run_locally, converter) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% processing gff-file no 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prediction_genes_fn = [predictions_dir '/anno/genes.mat'];
prediction_genes_dir = [predictions_dir '/anno/'];
fn_info = [predictions_dir '/anno/out'];
version = 'other';
seqid = '-';
source = '-';
level = 'level:gene:mRNA:CDS:five_prime_UTR,three_prime_UTR'
converter = '';

galaxy_gff2anno(genome_dir, predictions_fn, fn_info, prediction_genes_dir, version, seqid, source, level, ignore_missing_referenced_entries, run_locally, converter) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

label_format = 'Anno';
prediction_format =  'Anno';

fprintf('Performing evaluation\n') ;
galaxy_gene_eval(label_format, label_fn, label_genes_dir, prediction_format, predictions_fn, prediction_genes_dir, genome_dir, performance_fn, performance_dir,  eval_genomewide) ;

fprintf('Done.\n') ;


[ret, timedate] = unix('date') ;
fprintf('finished mgene_eval at %s', timedate) ;
