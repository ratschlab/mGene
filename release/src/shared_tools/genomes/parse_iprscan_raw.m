function genes = parse_iprscan_raw(iprscan_fname)
  %        genes = parse_iprscan_raw(iprscan_fname)


% count sequences

[tmp, num] = unix(sprintf('cut -f 1  %s | sort -u  |wc -l',iprscan_fname));
num_genes = str2double(num);
fprintf('%i sequences considered \n',num_genes);

% initiate gene structure

genes.name = '';
genes.crc64 = '';
genes.length = [];
genes.ipr_method = {};
genes.database_entry = {};
genes.database_desc = {};
genes.match_start = [];
genes.match_stop = [];
genes.evalue = [];
genes.status = [];
genes.IPR = {};
genes.IPR_desc = {};
genes.IPR_unique = {};
genes.GO ={};
genes.GO_desc = {};
genes.GO_unique ={};


genes(num_genes) = genes;

% start parsing
[fd msg] = fopen(iprscan_fname, 'r');
if fd<1, error('could not open file %s: %s',iprscan_fname, msg); end

genes_cnt = 0;
gene_names = {} ;
while ~feof(fd),
  elems = {} ;
  Line = fgetl(fd);
 
  if Line(1)=='#',
    continue ;
  end
 
  if ~ischar(Line),
    break ;
  end
  elems = separate(Line);
  assert(~isempty(elems{1}));
  
  identifier = elems{1};
  assert(~isempty(identifier));

  other_idx=strmatch(identifier,gene_names,'exact');
  if isempty(other_idx)
     genes_cnt =  genes_cnt +1;
     genes(genes_cnt).name = identifier;
     genes(genes_cnt).crc64 = elems{2};
     genes(genes_cnt).length = str2double(elems{3});;

     gene_names{end+1} = identifier ;
     gene_id = genes_cnt;
  else
    gene_id = other_idx;
  end
  assert(isequal(genes(gene_id).crc64,elems{2}));
  assert(isequal(genes(gene_id).length,str2double(elems{3})));

  genes(gene_id).ipr_method{end+1} = elems{4};
  genes(gene_id).database_entry{end+1} = elems{5};
  genes(gene_id).database_desc{end+1} = elems{6};
  genes(gene_id).match_start(end+1) = str2double(elems{7});
  genes(gene_id).match_stop(end+1) = str2double(elems{8});
  genes(gene_id).evalue{end+1} = elems{9};
  genes(gene_id).status{end+1} = elems{10};
  
  genes(gene_id).IPR{end+1} = elems{12};
  genes(gene_id).IPR_desc{end+1} = elems{13};
  
  if length(elems)>13
    genes(gene_id).GO_desc{end+1} = elems{14};
    GO_idx=strfind(elems{14},'(GO:');
    for i=GO_idx
      genes(gene_id).GO{end+1} = elems{14}(GO_idx+1:GO_idx+10);
    end
  else
    genes(gene_id).GO{end+1} = '';
    genes(gene_id).GO_desc{end+1} = '';
  end
  
  if mod(genes_cnt,1000)==0||genes_cnt==num_genes
    fprintf('  parsed %i genes\r', genes_cnt);  
  end
end 

% assert(genes_cnt==num_genes)
genes=genes(1:genes_cnt);

for i=1:length(genes)
  genes(i).IPR_unique=setdiff(unique(genes(i).IPR),'NULL');
  genes(i).GO_unique=setdiff(unique(genes(i).GO),'NULL');
end

fclose(fd);



