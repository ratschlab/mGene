function chunk_genes4iprscan(genes,genome_info, nb_per_chunk,fn)

nb_chunks = floor(length(genes)/nb_per_chunk);
idx = 0;
for i = 1:nb_chunks
  cds_fn = sprintf('%s_cds_chunk%i.fasta',fn,i);
  protein_fn = sprintf('%s_aa_chunk%i.fasta',fn,i);
  write_mrna_aa_fasta(genes(idx+1:idx+nb_per_chunk),genome_info,cds_fn,protein_fn)
  idx=idx+nb_per_chunk;
end
cds_fn = sprintf('%s_cds_chunk%i.fasta',fn,i+1);
protein_fn = sprintf('%s_aa_chunk%i.fasta',fn,i+1);
write_mrna_aa_fasta(genes(idx+1:end),genome_info,cds_fn,protein_fn)
