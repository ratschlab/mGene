function write_genes_cds2fasta(genes,genome_info,cds_fn,protein_fn, half_open_to_closed, append)

% write_fasta(genes,genome_info,cds_fn,protein_fn)
%   
% loads mRNA sequence for all genes and writes it to a fasta file (cds_fn).
% If protein_fn is specified, the sequences will also be translated and
% written to the fasta file protein_fn.
% 
% INPUT - genes required fields:  genes.name
%                                 genes.chr 
%                                 genes.strand
%                                 genes.cds_exons
%   
%        - genome_info
%        - cds_fn 
%        - protein_fn (optional)
%
% EXAMPLE 
% [a, b] = rgasp_genome_config_dir('japonica'); write_genes_cds2fasta(genes, init_genome(b), '~/tmp/jap_merged_gene_nuc', '~/tmp/jap_merged_gene_aa', 0)
%
% 
% SEE ALSO  write_fasta, load_genomic
num_arg = nargin;
if nargin==1
	paths
	PAR = genes;
	genes = PAR.genes;
	genome_info = PAR.genome_info;
	cds_fn = PAR.cds_fn;
	protein_fn = PAR.protein_fn;
	half_open_to_closed = PAR.half_open_to_closed;
	append = PAR.append;
	num_arg = 6;
end  

if num_arg>2
	if num_arg>=6 && append
		[fd1 msg] = fopen(cds_fn, 'a');
	else
		[fd1 msg] = fopen(cds_fn, 'w');
	end
else
  fd1 = 1;
end
if num_arg>3
  if num_arg>=6 && append
    [fd2 msg] = fopen(protein_fn, 'a');
  else
    [fd2 msg] = fopen(protein_fn, 'w');
  end
end
for i=1:length(genes)
  for j=1:length(genes(i).cds_exons)
    cds_exons = genes(i).cds_exons{j};
	if isempty(cds_exons)
		continue
	end
    if genes(i).strand=='+' && half_open_to_closed
      cds_exons(:,2) = cds_exons(:,2)-1;
    elseif half_open_to_closed
      cds_exons(:,1) = cds_exons(:,1)+1;
    end
    cds = load_genomic(genes(i).chr, char(genes(i).strand),cds_exons(:,1),cds_exons(:,2),genome_info,1);
    f = find_open_frames(cds,1) ;
	
	if ~(mod(length(cds),3)==0 && f(1)==1 && any(cds(1:3)=='atg')&&any(ismember({'taa', 'tag', 'tga'}, cds(end-2:end))))
		warning('found invalid ORF')
		mod(length(cds),3)
		f
		cds(1:3)
		cds(end-2:end)
		continue
	end
    %write_fasta(fd1, genes(i).name, cds) ;
    write_fasta(fd1, genes(i).transcripts{j}, cds) ;
    if exist('fd2')
      aa = translate(cds);
      %write_fasta(fd2, genes(i).name, aa) ;
      write_fasta(fd2, genes(i).transcripts{j}, aa) ;
    end
  end
end
fclose(fd1);
if num_arg>3
	fclose(fd2);
end
