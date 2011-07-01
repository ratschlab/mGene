function genes = parse_wetlab_txt(txt_fname)
% genes = parse_wetlab_txt(txt_fname)

[tmp,tmp2] =   unix(sprintf('grep GRseq %s| wc -l',txt_fname));

tmp2 = separate(tmp2,' ');
num_genes = str2num(tmp2{1});

genes = init_genes(num_genes);
[fd msg] = fopen(txt_fname, 'r');
genes_cnt = 0;
while ~feof(fd),
  
  Line = fgetl(fd);
  if strfind(Line,'reference = ')==1
    genes_cnt = genes_cnt+1;
    genes(genes_cnt).chr= Line(13:end);
    continue
  end
  elems = separate(Line);
  genes(genes_cnt).name = elems{2}(6:end); 
  exons = separate(elems{3},',');
  for j=1:length(exons)
    exon = separate(exons{j},'-');
    genes(genes_cnt).exons{1}(j,1) = str2num(exon{1});
    genes(genes_cnt).exons{1}(j,2) = str2num(exon{2});
  end
  genes(genes_cnt).start = min(min(genes(genes_cnt).exons{1}));
  genes(genes_cnt).stop = max(max(genes(genes_cnt).exons{1}));
end % end first parse
% assert(genes_cnt==num_genes)
genes=genes(1:genes_cnt);
fclose(fd);
fprintf('\n\n');
  
  