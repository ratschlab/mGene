function genes = parse_wetlab_psl(psl_fname,date)
% genes = parse_wetlab_txt(txt_fname)
  
if nargin<2
  date=[];
end
[tmp,tmp2] =   unix(sprintf('wc -l %s',psl_fname));

tmp2 = separate(tmp2,' ');
num_genes = str2num(tmp2{1});

genes = init_genes(num_genes);
[fd msg] = fopen(psl_fname, 'r');
genes_cnt = 0;
while ~feof(fd),
  
  Line = fgetl(fd);
  elems = separate(Line);
  genes_cnt = genes_cnt+1;
  genes(genes_cnt).date = date;
  genes(genes_cnt).chr = elems{14};
  genes(genes_cnt).strand = elems{9};
  genes(genes_cnt).name = elems{10}; 
  genes(genes_cnt).start = str2num(elems{16});
  genes(genes_cnt).stop = str2num(elems{17});
  genes(genes_cnt).match = str2num(elems{1});
  genes(genes_cnt).mismatch = str2num(elems{2});
  genes(genes_cnt).QgapCount = str2num(elems{5});
  genes(genes_cnt).QgapBases = str2num(elems{6});
  num_exons = str2num(elems{18});
  exons = separate(elems{21},',');
  exons_len = separate(elems{19},',');
  for j=1:num_exons 
    genes(genes_cnt).exons{1}(j,1) = str2num(exons{j});
    genes(genes_cnt).exons{1}(j,2) = str2num(exons{j})+str2num(exons_len{j});
  end
end % end first parse
% assert(genes_cnt==num_genes)
genes=genes(1:genes_cnt);
fclose(fd);
fprintf('\n\n');
 
  
  