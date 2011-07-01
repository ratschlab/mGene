function STATS = genes_statistics(genes,source)
% STATS = genes_statistics(genes)
  
if nargin<2
  source='mGene';
end
  
seqids = unique({genes.chr});  
STATS.gene = 0;
STATS.mRNA = 0;
STATS.CDS = 0;
STATS.five_prime_UTR = 0;
STATS.three_prime_UTR = 0;
for g=1:length(genes)
  gene=genes(g);
  
  STATS.gene = STATS.gene+1;
  STATS.mRNA  = STATS.mRNA+length(gene.transcripts);
  for j=1:length(gene.transcripts)
	if ~isempty(gene.cds_exons)
	    STATS.CDS = STATS.CDS+size(gene.cds_exons{j},1);
	    STATS.five_prime_UTR = STATS.five_prime_UTR +size(gene.utr5_exons{j},1);
	    STATS.three_prime_UTR = STATS.three_prime_UTR +size(gene.utr3_exons{j},1);
	end
  end
end

fprintf('\n------Genes SUMMARY-----\n');
fprintf('genes containes following seqids in field chr:\n');
for s=1:length(seqids)
  fprintf('%s\n',seqids{s})
end
fprintf('genes containes following sources:\n');
fprintf('%s\t%i\n',source,(STATS.gene+STATS.mRNA+STATS.five_prime_UTR+STATS.three_prime_UTR+STATS.CDS));
 
fprintf('\ngenes containes following types:\n');

types = fieldnames(STATS);
for t=1:length(types)
  fprintf('%s\t%i\n',types{t},STATS.(types{t}));
end
