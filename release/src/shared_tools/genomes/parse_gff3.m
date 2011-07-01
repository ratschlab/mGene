function [genes, FAULTS] = parse_gff3(gff_fname, save_fname, thesources, contig, TYPE,log_fname)
%        genes = parse_gff3(gff_fname, save_fname, thesources, contig, TYPE.log_fname)
% can also be used for gff2 format:
% TYPE.short_names = {'CDS','exon','coding_exon'};
% TYPE.IDs = {'1' '2' '3'};
% TYPE.hierachies = [2 3 3];
  
if nargin<3
  thesources = [];
end
  
  
if nargin<4
  contig = [];
end

if nargin<5||isempty(TYPE)
  %%% SOFA mapping
  SOFA_short_names = {'gene','mRNA','exon','CDS','intron','polyA_sequence','polyA_site','five_prime_UTR','three_prime_UTR'};
  SOFA_IDs = {'SO:0000704','SO:0000234','SO:0000147','SO:0000316','SO:0000188','SO:0000610','SO:0000553','SO:0000204','SO:0000205'};
  hierachies = [1 2 3 3 3 3 3 3 3];
else 
  SOFA_short_names = TYPE.short_names;
  SOFA_IDs = TYPE.IDs;
  hierachies = TYPE.hierachies;
end
 
status={'Predicted','Partially_confirmed','Confirmed'};

%%% GFF file summary

filtered_fname = gff3_statistics(gff_fname, thesources, contig);


if ~isempty(save_fname) && fexist(sprintf('%s.mat', save_fname))
  warning(sprintf('file %s.mat already exists\n', save_fname));
  %fprintf('\nwill be overwritten %s\n', save_fname);
  %load(save_fname,'genes');
  %assert(length(genes)==num_genes) 
end

FAULTS.Gene.no_strand = {};
FAULTS.Gene.no_id = {};
FAULTS.Gene.id_not_unique = {};
FAULTS.Transcript.unknown_parent = {};
FAULTS.Transcript.no_gene = {};
FAULTS.Transcript.diff_strand_contig = {};
FAULTS.Transcript.id_not_unique = {};
FAULTS.Transcript.outside_boundaries = {};
FAULTS.Exon.no_transcript = {};
FAULTS.Exon.diff_strand_contig = {};


% first parse: genes only are extracted 

types = SOFA_short_names(hierachies==1);
types_IDs = SOFA_IDs(hierachies==1);

fprintf('\n------starting first parse: only look for the following types\n');

tmp_fname1 = tempname;
unix(sprintf('touch %s',tmp_fname1));
for i=1:length(types)
  [tmp,num] = unix(sprintf('cut -f 3 %s | grep ^%s\t |wc -l',filtered_fname,types{i}));
  unix(sprintf('grep \t%s\t %s >> %s',types{i},filtered_fname,tmp_fname1));
  num_genes(i) = str2double(num);
  fprintf('%s\t%i\n',types{i},num_genes(i));
end

num_genes = sum(num_genes);
if num_genes ==0
  warning('no genes found; looking for hierachy 2 (mRNAs instead)\n')
  types = SOFA_short_names(hierachies==2);
  types_IDs = SOFA_IDs(hierachies==2);
  for i=1:length(types)
    [tmp,num] = unix(sprintf('cut -f 3 %s | grep ^%s\t |wc -l',filtered_fname,types{i}));
    unix(sprintf('grep \t%s\t %s >> %s',types{i},filtered_fname,tmp_fname1));
    num_genes(i) = str2double(num);
    fprintf('\t%s\t%i\n',types{i},num_genes(i));
  end
  num_genes
  num_genes = sum(num_genes);
end
if num_genes ==0
   warning('no genes found')
  return
end

genes = init_genes(num_genes);
[fd msg] = fopen(tmp_fname1, 'r');

if fd<1, error('could not open file %s: %s',tmp_fname1, msg); end

genes_cnt = 0;
gene_names = {} ;
while ~feof(fd),

  attributes = {} ;
  type_name = '' ;
  elems = {} ;
  
  
  Line = fgetl(fd);
  elems = separate(Line);
  assert(~isempty(elems{1}));
  
  if Line(1)=='#',
    continue ;
  end
  if ~isempty(contig)
    if ~isequal(upper(elems{1}),upper(contig)),
      continue ;
    end
  end
  if ~ischar(Line),
    break ;
  end
  assert(str2double(elems{4})<=str2double(elems{5}) || str2double(elems{5})<0)
  if~(elems{7}=='+' || elems{7}=='-') ;
    if elems{7}~='.',
      FAULTS.Gene.no_strand{end+1} = elems{9};
      continue
      % warning('strand missing')
    end ;
  end
  
  if ~isempty(thesources)
    if isempty(strmatch(elems{2}, thesources))
      error('something wrong with source')
    end  
  end
  
  %if ~ismember(elems{3},types) && ~ismember(elems{3},types_IDs)
  if isempty(strmatch(elems{3}, {types{:} types_IDs{:}}))
    continue ;
  end
  attributes = separate(elems{9}, ';');
  
  if ismember(elems{3},types);
    type_name = elems{3};
  elseif ismember(elems{3},types_IDs);
    type_name =types(isequal(elems{3},types_IDs));
  else
    error('wrong type')
  end
  

  identifier = '';    alias = '';
  for i=1:length(attributes),
    id=strfind(attributes{i},'ID=');
    if ~isempty(id),
      identifier = attributes{i}(id+3:end);
      break;
      % keyboard
    end
  end
  if isempty(identifier)
    FAULTS.Gene.no_id{end+1} = elems{9};
    % warning('no ID found 1')
    identifier = attributes{1};
  end
      
  for i=1:length(attributes) 
    id=strfind(attributes{i},'Alias=');
    if ~isempty(id),
      alias = attributes{i}(id+6:end);
      break;
    end 
  end
  if length(identifier)>=3 
    assert(~isequal(identifier(1:3),'ID=')) ;
  end
  assert(~isempty(identifier));

  other_idx=strmatch(identifier,gene_names,'exact');
  if ~isempty(other_idx)
    % warning('gene appears more than once\n')
    FAULTS.Gene.id_not_unique{end+1} = elems{9}; 
    assert(isequal(genes(other_idx).alias,alias));
    assert(isequal(genes(other_idx).chr,elems{1}));
    assert(isequal(genes(other_idx).strand,elems{7})) ;
    if ~isequal(genes(other_idx).start,str2double(elems{4}))
      genes(other_idx).start = min(genes(other_idx).start,str2double(elems{4}));
    end
    if ~isequal(genes(other_idx).stop,str2double(elems{5}))
      genes(other_idx).stop = max(genes(other_idx).stop,str2double(elems{5}));
    end 
    % assert(isequal(genes(other_idx).start,str2double(elems{4})));
    % assert(isequal(genes(other_idx).stop,str2double(elems{5})));
    continue
  end
  genes_cnt =  genes_cnt +1;
  genes(genes_cnt).name = identifier;
  genes(genes_cnt).alias = alias;
  genes(genes_cnt).chr = elems{1};
  genes(genes_cnt).strand = elems{7} ;
  genes(genes_cnt).start(end+1) = str2double(elems{4});
  genes(genes_cnt).stop(end+1) = str2double(elems{5});
  genes(genes_cnt).is_correctly_gff3_referenced = 1 ;

  gene_names{end+1} = identifier ;
  if mod(genes_cnt,1000)==0||genes_cnt==num_genes
    fprintf('  parsed %i genes\r', genes_cnt);  
  end
end % end first parse
% assert(genes_cnt==num_genes)
genes=genes(1:genes_cnt);
fclose(fd);

fprintf('\n\n');
if length(FAULTS.Gene.no_strand)>0 || length(FAULTS.Gene.no_id)>0 || length(FAULTS.Gene.id_not_unique)>0,
  fprintf('Following problems were found:\n') ;
else
  fprintf('No problems found.\n') ;
end ;
if length(FAULTS.Gene.no_strand)>0,
  fprintf('number of genes without annotated strand (column 7):\t\t %i\n', length(FAULTS.Gene.no_strand)) ;
  fprintf('\t=> genes are ignored\n') ;
end 
if length(FAULTS.Gene.no_id)>0,
  fprintf('number of genes without annotated ID (column 9, ID=... ):\t %i\n', length(FAULTS.Gene.no_id)) ;
  fprintf('\t=> first part (to the first semicolon) of attribute field (column 9) is used as ID\n') ;
end
if length(FAULTS.Gene.id_not_unique)>0,
  fprintf('number of genes with non-unique ID (column 9, ID=... ):\t\t %i\n', length(FAULTS.Gene.id_not_unique)) ;
  fprintf('\t=> all on the same seqid and strand, min of start (column 4) and max of end (column 5) used\n ') ;
end ;

unix(sprintf('rm %s',tmp_fname1));




% second parse: mRNAs only are extracted 

types = SOFA_short_names(hierachies==2);
types_IDs = SOFA_IDs(hierachies==2);

fprintf('\n------starting second parse: only look for the following types\n');
tmp_fname1 = tempname;
unix(sprintf('touch %s',tmp_fname1));
num_transcripts = zeros(1,length(types)) ;
for i=1:length(types)
  [tmp,num] = unix(sprintf('cut -f 3 %s | grep %s\t |wc -l',filtered_fname,types{i}));
  unix(sprintf('grep \t%s\t %s >> %s',types{i},filtered_fname,tmp_fname1));
  num_transcripts(i) = str2double(num);
  fprintf('%s\t%i\n',types{i},num_transcripts(i));
end

num_transcripts = sum(num_transcripts);
transcript_gene_map = zeros(1,num_transcripts);
transcript_names = cell(1,num_transcripts);
transcript_ids = cell(1,num_transcripts);
gene_names = {genes.name} ;

[fd msg] = fopen(tmp_fname1, 'r');
transcript_cnt = 0;
while ~feof(fd),

  attributes = {} ;
  type_name = '' ;
  elems = {} ;
  
  
  Line = fgetl(fd);
  elems = separate(Line);
  assert(~isempty(elems{1}));
  
  if Line(1)=='#',
    continue ;
  end
  if ~isempty(contig)
    if ~isequal(upper(elems{1}),upper(contig)),
      continue ;
    end
  end
  if ~ischar(Line),
    break ;
  end
  assert(str2double(elems{4})<=str2double(elems{5}) || str2double(elems{5})<0)
  if~(elems{7}=='+' || elems{7}=='-') ;
    if elems{7}~='.',
      % warning('strand missing')
    end ;
  end
  
  if ~isempty(thesources)
    if isempty(strmatch(elems{2}, thesources))
      error('something wrong with source')
    end  
  end
  
  if ~ismember(elems{3},types) && ~ismember(elems{3},types_IDs)
    continue ;
  end
  attributes = separate(elems{9}, ';');
  
  if ismember(elems{3},types);
    type_name = elems{3};
  elseif ismember(elems{3},types_IDs);
    type_name =types(isequal(elems{3},types_IDs));
  else
    error('wrong type')
  end

  identifier = ''; parent = ''; pred_status = '';
  for i=1:length(attributes),
    id=strfind(attributes{i},'ID=');
    if ~isempty(id),
      identifier = attributes{i}(id+3:end);
    end 
    id=strfind(upper(attributes{i}),'PARENT=');
    if ~isempty(id),
      parents=separate(attributes{i}(id+7:end), ',') ;
      parent =parents{1} ;
    end
    id=strfind(attributes{i},'prediction_status=');
    if ~isempty(id),
      pred_status = attributes{i}(id+18:end);
      pred_status = strmatch(pred_status,status)-2;
    end
  end
  if isempty(identifier)&&isempty(parent)% &length(attributes)==1
    identifier = attributes{1};
    parent = attributes{1};
  end    
  parent_idx=strmatch(parent,gene_names,'exact');
  if isempty(parent_idx)
    parent_idx=strmatch(identifier,gene_names,'exact');
  end
  if isempty(parent_idx),
    if isempty(parent),
      FAULTS.Transcript.unknown_parent{end+1} = elems{9};
      % warning(sprintf('Empty parent: using transcript name as parent (%s)', identifier)) ;
      parent = identifier ;
    end 
    %keyboard
    % warning(sprintf('Could not find parent gene ("%s"): Created a new gene descriptor', parent)) ;
    FAULTS.Transcript.no_gene{end+1} = elems{9};
    genes(end+1).name = parent ;
    genes(end).alias = parent ;
    genes(end).chr = elems{1};
    genes(end).strand = elems{7} ;
    genes(end).start(end+1) = 0 ;
    genes(end).stop(end+1) = inf ;
    genes(end).is_correctly_gff3_referenced = 0 ;
    gene_names = {genes.name} ;
    parent_idx=length(genes) ;
  end ;
  if ~isequal(elems{1}, genes(parent_idx).chr) ;
    FAULTS.Transcript.diff_strand_contig{end+1} = elems{9};
  end
  if ~isequal(elems{7}, genes(parent_idx).strand) ;
    FAULTS.Transcript.diff_strand_contig{end+1} = elems{9};
    % warning(sprintf('strand of parent gene (%s) different',genes(parent_idx).name))
    continue
  end
  assert(~isempty(parent_idx))
  if ~isempty(genes(parent_idx).transcripts)
    other_idx = strmatch(identifier,genes(parent_idx).transcripts,'exact');
    if ~isempty(other_idx);
      FAULTS.Transcript.id_not_unique{end+1} = elems{9};
      % warning('transcript appears more than once!')
      if ~isequal(genes(parent_idx).start(other_idx), str2double(elems{4})) 
        genes(parent_idx).start = min(genes(parent_idx).start(other_idx),str2double(elems{4}));             
      end 
      if ~isequal(genes(parent_idx).stop(other_idx), str2double(elems{5}))            
        genes(parent_idx).stop(other_idx) = max(genes(parent_idx).stop(other_idx),str2double(elems{5}));             
      end
      continue
    end
  end
  genes(parent_idx).transcripts{end+1} = identifier;
  if ~isempty(pred_status)
    genes(parent_idx).transcript_status(end+1) = pred_status;
  end
  if isempty(genes(parent_idx).exons) 
    genes(parent_idx).exons = cell(1,1) ;
  else
    genes(parent_idx).exons{end+1} = [] ;
  end
  %transcript must lie within gene boundaries
  if ~(str2double(elems{4})>=genes(parent_idx).start(1)&&str2double(elems{4})<=genes(parent_idx).stop(1))
    FAULTS.Transcript.outside_boundaries{end+1} = elems{9};
    % warning(sprintf('transcript outside gene boundaries (%s; %s)', identifier, parent)) ;
  end
  if ~(str2double(elems{5})>=genes(parent_idx).start(1)&&str2double(elems{5})<=genes(parent_idx).stop(1))
    FAULTS.Transcript.outside_boundaries{end+1} = elems{9};
    % warning(sprintf('transcript outside gene boundaries (%s; %s)', identifier, parent)) ;
  end
  genes(parent_idx).start(end+1) = str2double(elems{4});
  genes(parent_idx).stop(end+1) = str2double(elems{5}); 
  transcript_cnt = transcript_cnt +1;
  transcript_gene_map(transcript_cnt) = parent_idx;
  transcript_names{transcript_cnt} = identifier;
  tmp = separate(identifier,':');
  if length(tmp)>1
    transcript_ids{transcript_cnt} = tmp{2}; 
  else
    transcript_ids{transcript_cnt} = tmp{1}; 
  end
  if mod(transcript_cnt,1000)==0||transcript_cnt == num_transcripts
    fprintf('  parsed %i transcripts\r', transcript_cnt);  
  end
end %end second parse
transcript_names = transcript_names(1:transcript_cnt); 
transcript_ids = transcript_ids(1:transcript_cnt); 

% assert(transcript_cnt==num_transcripts)
fclose(fd);
unix(sprintf('rm %s',tmp_fname1));

fprintf('\n\n');
if length(FAULTS.Transcript.unknown_parent)>0 || length(FAULTS.Transcript.no_gene)>0 || length(FAULTS.Transcript.diff_strand_contig)>0 || length(FAULTS.Transcript.id_not_unique)>0 || length(FAULTS.Transcript.outside_boundaries)
  fprintf('Following problems were found:\n')
else
  fprintf('No problems found.\n')
end ;
if length(FAULTS.Transcript.unknown_parent)>0 
  fprintf('Number of transcripts without annotated parent (keyword ''Parent='' missing in column 9):\t %i\n',length(FAULTS.Transcript.unknown_parent))
  fprintf('\t=> using ID as Parent\n')
end
if length(FAULTS.Transcript.no_gene)>0
  fprintf('Number of transcripts without annotated parent gene:\t\t\t %i\n',length(FAULTS.Transcript.no_gene))
  fprintf('\t=> created a new gene descriptor\n')
end
if length(FAULTS.Transcript.diff_strand_contig)>0
  fprintf('Number of transcripts with different strand or contig than parent gene:\t %i\n',length(FAULTS.Transcript.diff_strand_contig))
  fprintf('\t=> transcript ignored\n')
end
if length(FAULTS.Transcript.id_not_unique)>0
  fprintf('Number of transcripts with non-unique ID (column 9, ID=... ):\t\t %i\n',length(FAULTS.Transcript.id_not_unique))
  fprintf('\t=> only one transcript added, min of start (column 4) and max of end (column 5) used\n')
end
if length(FAULTS.Transcript.outside_boundaries)
  fprintf('Number of transcripts outside boundaries of gene:\t\t\t %i\n',length(FAULTS.Transcript.outside_boundaries))
  fprintf('\t=> gene boundaries extended\n')
end ;

% correcting gene boundaries for orphan transcripts
for i=genes_cnt+1:length(genes),
  genes(i).start(1) = min(genes(i).start(2:end)) ;
  genes(i).stop(1) = max(genes(i).stop(2:end)) ;
end ;

% third parse: genes are filled with exons

% keyboard
types = SOFA_short_names(hierachies==3);
types_IDs = SOFA_IDs(hierachies==3);
fprintf('\n------starting third parse: only look for the following type\n');
tmp_fname1 = tempname;
unix(sprintf('touch %s',tmp_fname1));

for i=1:length(types)
  [tmp,num] = unix(sprintf('cut -f 3 %s | grep %s\t |wc -l',filtered_fname,types{i}));
  unix(sprintf('grep \t%s\t %s >> %s',types{i},filtered_fname,tmp_fname1));
  num_segments(i) = str2num(num);
  fprintf('%s\t%i\n',types{i},num_segments(i));
end
num_segments = sum(num_segments);

%num_segments = zeros(1,length(types)) ;
%for i=1:length(types)
%  [tmp,num] = unix(sprintf('cut -f 3 %s | grep %s |wc -l',tmp_fname,types{i}));
%  num_segments(i) = str2double(num) ;
%  fprintf('%s\t%i\n',types{i},num_segments(i));
%end
%num_segments = sum(num_segments);

%disp('creating hash table for transcript names') ;
%[transcript_hash_names, transcript_hash_idx]= init_hash_search(transcript_names) ;

[fd msg] = fopen(tmp_fname1, 'r');
disp(msg);
segment_cnt = 0;
segment_skipped=0;
while ~feof(fd),

  attributes = {} ;
  type_name = '' ;
  elems = {} ;
  
  Line = fgetl(fd);
  elems = separate(Line);
  assert(~isempty(elems{1}));
  
  if Line(1)=='#',
    continue ;
  end
  if ~isempty(contig)
    if ~isequal(upper(elems{1}),upper(contig)),
      continue ;
    end
  end
  if ~ischar(Line),
    break ;
  end
  assert(str2double(elems{4})<=str2double(elems{5}) || str2double(elems{5})<0)
  if~(elems{7}=='+' || elems{7}=='-') ;
    if elems{7}~='.',
      % warning('strand missing')
    end ;
  end
  
  if ~isempty(thesources)
    if isempty(strmatch(elems{2}, thesources))
      error('something wrong with source')
    end  
  end
  
  if ~ismember(elems{3},types) && ~ismember(elems{3},types_IDs)
    continue ;
  end
  attributes = separate(elems{9}, ';');
  
  if ismember(elems{3},types);
    type_name = elems{3};
  elseif ismember(elems{3},types_IDs);
    type_name =types(isequal(elems{3},types_IDs));
  else
    error('wrong type')
  end

  switch type_name ,
   case {'CDS','coding_exon'}
    type = 3;
   case 'three_prime_UTR'
    type = 4;
   case 'five_prime_UTR'
    type = 2; 
   case 'exon'
    type = 0;
   otherwise
    segment_skipped = segment_skipped+1;
    % disp(Line);
    continue
  end
  
  parent = '';
  for i=1:length(attributes), 
    id=strfind(upper(attributes{i}),'PARENT=');
    if ~isempty(id),
      parent=attributes{i}(id+7:end) ;
      %parent = parents{1} ;
    end
  end
  if isempty(parent) %&  length(attributes)==1
    parent = attributes{1};
  end
  assert(~isempty(parent));
  parents = separate(parent, ',') ;
  for parent=parents,

    if ~ischar(parent),
      assert(length(parent)==1) ;
      parent=parent{1} ;
    end ;

    parent_idx = strmatch(parent, transcript_names, 'exact');
    %parent_idx = hash_search(transcript_hash_names, transcript_hash_idx, parent) ;
    %assert(isequal(parent_idx, parent_idx_)) ;

    id_seg = separate(parent, ':') ;
    if length(id_seg)>1
      id_seg = id_seg{2} ;
    end
    if isempty(parent_idx)
      id_seg = separate(parent, ':') ;
      if length(id_seg)>1
        id_seg = id_seg{2} ;
      end
      parent_idx=strmatch(id_seg, transcript_ids, 'exact');
    end
    gene_idx = [];
    if isempty(parent_idx),
      gene_idx =strmatch(parent,gene_names,'exact');
    end

    if isempty(parent_idx)&&isempty(gene_idx),
      % make sure that we only create a parent, if none of the parents exists
      assert(~isempty(parent)) ;
      % warning(sprintf('Could not find parent transcript ("%s"): Created a new gene and transcript descriptor', parent)) ;
      FAULTS.Exon.no_transcript{end+1} = elems{9};
      genes(end+1).name = [parent '-gene'] ;
      genes(end).alias = parent ;
      genes(end).chr = elems{1};
      genes(end).strand = elems{7} ;
      genes(end).start(end+1) = 0 ;
      genes(end).stop(end+1) = inf ;
      genes(end).is_correctly_gff3_referenced = 0 ;
      genes(end).transcripts{end+1} = parent ;
      genes(end).exons{end+1} = [] ;
      
      transcript_cnt = transcript_cnt +1;
      transcript_gene_map(transcript_cnt) = length(genes) ;
      transcript_names{transcript_cnt} = parent ;

      % add to hash table
      %[transcript_hash_names, transcript_hash_idx]= init_hash_search(transcript_names, transcript_cnt, transcript_hash_names, transcript_hash_idx) ;
      break ;
    else
      break ;
    end ;
  end ;

  for parent = parents,
    segment_invalid = 0 ;
    if ~ischar(parent),
      parent=parent{1} ;
    end ;

    parent_idx = strmatch(parent, transcript_names, 'exact') ;
    %parent_idx = hash_search(transcript_hash_names, transcript_hash_idx, parent) ;
    %assert(parent_idx==parent_idx_) ;

    if isempty(parent_idx),
      fprintf('parent %s not found\n', parent) ;
      gene_idx = transcript_gene_map(parent_idx);
      if isempty(gene_idx),
        continue ;
      end
      genes(gene_idx).is_correctly_gff3_referenced = 0 ;
      transcript_idx = 1;
    else
      assert(~isempty(parent_idx) && length(parent_idx)==1);
      gene_idx = transcript_gene_map(parent_idx);
      transcript_idx = strmatch(parent, genes(gene_idx).transcripts, 'exact');
    end ;
    
    if isempty(transcript_idx)
      id_seg = separate(parent, ':') ;
      if length(id_seg)>1
        id_seg = id_seg{2} ;
      end
      for j = 1:length(genes(gene_idx).transcripts)
        tmp =  separate(genes(gene_idx).transcripts{j}, ':') ;
        id_trans{j} = tmp{2};
      end
      transcript_idx=strmatch(id_seg,id_trans,'exact');
    end
    assert(~isempty(transcript_idx) & length(transcript_idx)==1);

    if ~isequal(elems{1}, genes(gene_idx).chr)
      segment_invalid = 1 ;
      FAULTS.Exon.diff_strand_contig{end+1} = elems{9};
      % fprintf('Segment''s contig name does not match parent (is %s, should be %s). Ignoring segment.\nLine:\n %s\n', elems{1}, genes(gene_idx).chr, Line) ;
    end ;
    if ~isequal(elems{7}, genes(gene_idx).strand),
      segment_invalid = 1 ;
      FAULTS.Exon.diff_strand_contig{end+1} = elems{9};
      % fprintf('Segment''s strand name does not match parent (is %s, should be %s). Ignoring segment.\nLine:\n %s\n', elems{7}, genes(gene_idx).strand, Line) ;
    end ;

    if ~segment_invalid,
      if elems{7}=='+'
        genes(gene_idx).exons{transcript_idx}(end+1,:) = [str2double(elems{4}), str2double(elems{5})+1, type];
      else
        genes(gene_idx).exons{transcript_idx}(end+1,:) = [str2double(elems{4})-1, str2double(elems{5}), type];
      end
    end ;
%    assert(str2double(elems{4})>=genes(gene_idx).start(1+transcript_idx)&str2double(elems{4})<=genes(gene_idx).stop(1+transcript_idx))
%    assert(str2double(elems{5})>=genes(gene_idx).start(1+transcript_idx)&str2double(elems{5})<=genes(gene_idx).stop(1+transcript_idx))
  end ;
  segment_cnt = segment_cnt+1;
  if mod(segment_cnt,5000)==0
    fprintf('  parsed %i gene segments\r', segment_cnt);
  end
end % end third parse
fclose(fd);


fprintf('\n\n');
if length(FAULTS.Exon.no_transcript)>0 || length(FAULTS.Exon.diff_strand_contig)>0,
  fprintf('Following problems were found:\n');
else
  fprintf('No problems found.\n');
end ;
if length(FAULTS.Exon.no_transcript)>0, 
  fprintf('number of gene segments without parent transcript:\t %i\n', length(FAULTS.Exon.no_transcript)) ;
  fprintf('\t=> created a new gene and transcript descriptor\n') ;
end ;
if length(FAULTS.Exon.diff_strand_contig)>0,
  fprintf('number of gene segments with different strand or contig than parent transcript:\t %i\n', length(FAULTS.Exon.diff_strand_contig)) ;
  fprintf('\t=> Segment ignored\n') ;
end ;

for i=1:length(genes)
  genes(i).start = min(genes(i).start);
  genes(i).stop = max(genes(i).stop);
end 

if ~isequal(gff_fname, tmp_fname1),
  unix(['rm ' tmp_fname1]);
end ;

if ~isempty(save_fname),
   save(save_fname, 'genes', 'FAULTS');
end ;

