function [gene_count,genecover,exon_count,intron_count,transcript_count] = count_genestats(genesfile)
%function [gene_count,genecover,exon_count,intron_count,transcript_count] = count_genestats(genesfile)
%
% count:
% - the number of genes
% - the number of nucleotides covered by genes
% - the number of exons and introns in the splice graph
% - the number of transcripts

load(genesfile,'genes');

gene_count = length(genes);
genecover = 0;
exon_count = 0;
intron_count = 0;
transcript_count = 0;

for ix = 1:gene_count
  if mod(ix,1000)==0, fprintf('.');end;
  genelen = genes(ix).stop - genes(ix).start;
  assert(genelen>0)
  genecover = genecover + genelen;
  exon_count = exon_count + size(genes(ix).splicegraph{1},2);
  intron_count = intron_count + sum(sum(triu(genes(ix).splicegraph{2})));
  transcript_count = transcript_count + length(genes(ix).transcripts);
end

fprintf('\n');
