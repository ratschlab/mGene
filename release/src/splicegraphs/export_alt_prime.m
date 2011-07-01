disp('Loading genes');
load Data/elegans.genes.graph.mat
buffer_length = 0;

%fsock = 1;

fsock = fopen('Data/export_alt_prime.txt','w');

disp('Exporting 5 prime');
fprintf(fsock,'Alternative 5 prime exon locations\n');

for ix = 1:length(idx_alt_5prime)
  start = genes(idx_alt_5prime(ix)).splicegraph{1}(1, ...
						   exon_alt_5prime(ix));
  stop = genes(idx_alt_5prime(ix)).splicegraph{1}(2, ...
						   exon_alt_5prime(ix));
  chr = genes(idx_alt_5prime(ix)).chr;
  strand = genes(idx_alt_5prime(ix)).strands(1);
  name = genes(idx_alt_5prime(ix)).name;
  take_idx = find(name=='-');
  name = name(take_idx(1)+1:take_idx(2)-1);
  
  fprintf(fsock,'%s\t%s\t%c\t%d\t%d\n',name,chr,strand,...
	  start-buffer_length,stop+buffer_length);
end




disp('Exporting 3 prime');
fprintf(fsock,'\n\n\nAlternative 3 prime exon locations\n');

for ix = 1:length(idx_alt_3prime)
  start = genes(idx_alt_3prime(ix)).splicegraph{1}(1, ...
						   exon_alt_3prime(ix));
  stop = genes(idx_alt_3prime(ix)).splicegraph{1}(2, ...
						   exon_alt_3prime(ix));
  chr = genes(idx_alt_3prime(ix)).chr;
  strand = genes(idx_alt_3prime(ix)).strands(1);
  name = genes(idx_alt_3prime(ix)).name;
  take_idx = find(name=='-');
  name = name(take_idx(1)+1:take_idx(2)-1);

  fprintf(fsock,'%s\t%s\t%c\t%d\t%d\n',name,chr,strand,...
	  start-buffer_length,stop+buffer_length);
end

fclose(fsock);

