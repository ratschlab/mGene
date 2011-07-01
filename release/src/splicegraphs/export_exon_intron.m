%disp('Loading genes');
%load ../Data/elegans.genes.graph.mat
buffer_length = 0

%fsock = 1;

fsock = fopen('Data/export_exon_intron.txt','w');

disp('Exporting Exons in Intronic regions');
fprintf(fsock,'Alternative exon locations\n');

for ix = 1:length(idx_exon_intron)
  start = genes(idx_exon_intron(ix)).splicegraph{1}(1, ...
						   exon_exon_intron(ix));
  stop = genes(idx_exon_intron(ix)).splicegraph{1}(2, ...
						   exon_exon_intron(ix));
  chr = genes(idx_exon_intron(ix)).chr;
  strand = genes(idx_exon_intron(ix)).strands(1);
  name = genes(idx_exon_intron(ix)).name;
  take_idx = find(name=='-');
  name = name(take_idx(1)+1:take_idx(2)-1);
  
  fprintf(fsock,'%s\t%s\t%c\t%d\t%d\n',name,chr,strand,...
	  start-buffer_length,stop+buffer_length);
end




fclose(fsock);

