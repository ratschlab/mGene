function genes = parse_gff3_pseudogenes(gff_fname, save_fname,source,contig)
% Genes = parse_gff3(gff_fname, save_fname,source)
  
  
if nargin<3
  source = [];
end                                                                                     
if nargin<4                                                                           
  contig = [];                                                                            
end

%%% SOFA mapping
SOFA_short_names = {'pseudogene','exon','intron'};
SOFA_IDs = {'SO:0000704','SO:0000234','SO:0000147','SO:0000316','SO:0000188','SO:0000610','SO:0000553','SO:0000204','SO:0000205'};
hierachies = [1 3 3];
status={'Predicted','Partially_confirmed','Confirmed'};

%%% GFF file summary

tmp_fname = gff3_statistics(gff_fname, source,contig);

[tmp,num_genes] = unix(sprintf('cut -f 3 %s | grep gene |wc -l',tmp_fname));
num_genes = str2num(num_genes);

if fexist(sprintf('%s.mat',save_fname))
  fprintf('file %s already exists\n',save_fname)
  fprintf('\nloaded gene data from %s\n', save_fname);
  load(save_fname,'genes');
  %assert(length(genes)==num_genes) 
  return
end

% first parse: genes only are extracted 

types = SOFA_short_names(hierachies==1);
types_IDs = SOFA_IDs(hierachies==1);
fprintf('\n------starting first parse: only look for types\n');
for i=1:length(types)
  [tmp,num] = unix(sprintf('cut -f 3 %s | grep %s |wc -l',tmp_fname,types{i}));
  num_genes(i) = str2num(num);
  fprintf('%s\t%i\n',types{i},num_genes(i));
end
num_genes = sum(num_genes);
if num_genes ==0
  warning('no genes found; looking for hierachy 2 (mRNAs instead)\n')
  types = SOFA_short_names(hierachies==2);
  types_IDs = SOFA_IDs(hierachies==2);
  fprintf('\n------starting first parse: only look for types\n');
  for i=1:length(types)
    [tmp,num] = unix(sprintf('cut -f 3 %s | grep %s |wc -l',tmp_fname,types{i}));
    num_genes(i) = str2num(num);
    fprintf('%s\t%i\n',types{i},num_genes(i));
  end
  num_genes = sum(num_genes);
end
if num_genes ==0
   warning('no genes found')
  return
end
genes = init_genes(num_genes);
[fd msg] = fopen(tmp_fname, 'r');
genes_cnt = 0;
while ~feof(fd),
  parse_starter
  identifier = '';    alias = '';
  for i=1:length(attributes),
    id=strfind(attributes{i},'ID=');
    if ~isempty(id),
      identifier = attributes{i}(id+3:end);
    end 
  end
  assert(~isempty(identifier)) ;
  other_idx=strmatch(identifier,{genes(1:genes_cnt).name},'exact');
  if ~isempty(other_idx)
    warning('gene appears more than once\n')
    assert(isequal(genes(other_idx).chr,elems{1}));
    assert(isequal(genes(other_idx).strand,elems{7})) ;
    %assert(isequal(genes(other_idx).start,str2num(elems{4})));
    %assert(isequal(genes(other_idx).stop,str2num(elems{5})));
    genes(other_idx).start=min(genes(other_idx).start,str2num(elems{4})) ;
    genes(other_idx).stop=max(genes(other_idx).stop,str2num(elems{5})) ;
    continue
  end
  genes_cnt =  genes_cnt +1;
  genes(genes_cnt).name = identifier;
  genes(genes_cnt).transcripts{1} = identifier;  
  genes(genes_cnt).chr = elems{1};
  genes(genes_cnt).strand = elems{7} ;
  genes(genes_cnt).start(end+1) = str2num(elems{4});
  genes(genes_cnt).stop(end+1) = str2num(elems{5});
  genes(genes_cnt).exons = cell(1,1) ;
  fprintf('  parsed %i genes\r', genes_cnt);  
end % end first parse
% assert(genes_cnt==num_genes)
genes=genes(1:genes_cnt);
fclose(fd);
fprintf('\n\n');
gene_names = {genes.name} ;

% third parse: genes are filled with exons

types = SOFA_short_names(hierachies==3);
types_IDs = SOFA_IDs(hierachies==3);
fprintf('\n------starting third parse: only look for type\n');
for i=1:length(types)
  [tmp,num] = unix(sprintf('cut -f 3 %s | grep %s |wc -l',tmp_fname,types{i}));
  num_segments(i) = str2num(num);
  fprintf('%s\t%i\n',types{i},num_segments(i));
end
num_segments = sum(num_segments);

fprintf('skipping introns (redundant with exons)\n');

[fd msg] = fopen(tmp_fname, 'r');
disp(msg);
segment_cnt = 0;segment_skipped=0;

while ~feof(fd),

  parse_starter
  switch type_name ,
   case 'exon'
    type = 0;
   otherwise
    segment_skipped = segment_skipped+1;
    % disp(Line);
    continue
  end
  
  parent = '';
  for i=1:length(attributes), 
    id=strfind(attributes{i},'Name=');
    if ~isempty(id),
      parent = attributes{i}(id+5:end);
    end
  end
  assert(~isempty(parent));
  parents = separate(parent, ',') ;
  for parent=parents,
    parent_idx=strmatch(parent,gene_names,'exact');
    assert(~isempty(parent_idx) & length(parent_idx)==1);
    gene_idx = parent_idx;

    transcript_idx = strmatch(parent,genes(gene_idx).transcripts,'exact');
    assert(~isempty(transcript_idx) & length(transcript_idx)==1);

    assert(isequal(elems{1}, genes(gene_idx).chr)) ;
    assert(isequal(elems{7}, genes(gene_idx).strand)) ;

    if elems{7}=='+'
      genes(gene_idx).exons{transcript_idx}(end+1,:) = [str2num(elems{4}), str2num(elems{5})+1, type];
    else
      genes(gene_idx).exons{transcript_idx}(end+1,:) = [str2num(elems{4})-1, str2num(elems{5}), type];
    end
    assert(str2num(elems{4})>=genes(gene_idx).start&str2num(elems{4})<=genes(gene_idx).stop)
    assert(str2num(elems{5})>=genes(gene_idx).start&str2num(elems{5})<=genes(gene_idx).stop)
  end ;
  segment_cnt = segment_cnt+1;
  if mod(segment_cnt,100)==0
    fprintf('  parsed %i gene segments\r', segment_cnt);
  end
end % end third parse
fclose(fd);

assert(segment_cnt+segment_skipped==num_segments)
fprintf('\n\n');
for i=1:length(genes)
  genes(i).start = min(genes(i).start);
  genes(i).stop = max(genes(i).stop);
end

save(save_fname, 'genes');


