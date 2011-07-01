function genes = parse_gff2(gff_fname, save_fname, genome_info, source, contig, TYPE)
% genes = parse_gff2(gff_fname, save_fname, source, contig, type)
%
% parses gff2(gtf) files into the structure genes.
%
%   <gff_fname> is the full path to the gff(2) file to be parsed.
%
%   <save_fname> is the full path to the outfile to which the 
%    parsed data will be written. It doesn't have to be existent,
%    in case it will be overwritten.
%
%   <genome_info> is the genome_info structure fitting to the gff data
%
%   <source> is an optional parameter which contains the names of
%    accepted data-sources. All lines in which another source is
%    specified will be ignored. If not set, all lines will be
%    accepted. e.g. source={'source_1', 'source_2', ...}
%
%   <contig> is an optional parameter which should contain the 
%    names of the wanted FPC-IDs (or chromosome name). It is used
%    to filter the lines for wanted contigs. If not set, nodoc 
%    lines will be ignored. e.g. contig={'chr1', 'chr2', ...}
%   
%   <TYPE> is an optional parameter which is used to filter the
%    lines for wanted feature types (CDS, exon...). It has to be 
%    a structure with fields 
%    TYPE.short_names={'name_1', 'name_2', ...}.
%    If not set, no lines will be filtered out.



if nargin < 2
  warning('to few parameters');
end

if nargin < 3
  source=[];
end

if nargin < 4
  contig=[];
end

% no TYPEs specified
if nargin<5||isempty(TYPE)
  %%% SOFA mapping
  SOFA_short_names = {'exon','CDS','intron','polyA_sequence','polyA_site','five_prime_UTR','three_prime_UTR'};
  %SOFA_IDs = {'SO:0000147','SO:0000316','SO:0000188','SO:0000610','SO:0000553','SO:0000204','SO:0000205'};
  %hierachies = [3 3 3 3 3 3 3];
else 
  SOFA_short_names = TYPE.short_names;
  %SOFA_IDs = TYPE.IDs;
  %hierachies = TYPE.hierachies;
end

% test if infile is existent
if isempty(gff_fname)
  warning('no infile specified');
  return;
elseif ~fexist(sprintf('%s',gff_fname))
  warning('file not existent or maybe wrong path');
  return;
end

% test if file 'save_fname' already exists
if isempty(save_fname)
  warning('no outfile specified');
  return;
elseif fexist(sprintf('%s.mat',save_fname))
  warning(sprintf('file %s.mat already exists. Will be overwritten!\n',save_fname))
end

% test if gene_id is in correct column
[tmp,label] = unix(sprintf('cut -f 9 %s | cut -d \\" -f 1 | sort -u',gff_fname));
if ~strcmp(char('gene_id'), strtrim(label))
  warning(sprintf('file %s wrong formated\n',gff_fname))
  return
end

% get the number of genes in gff_file
% not yet filtered so num is to high if filters exist
[tmp,num] = unix(sprintf('cut -f 9 %s | cut -d \\" -f 2 | sort -u | wc -w',gff_fname));
num = int32(str2double(num));
fprintf('Unfiltered number of genes in infile: %u\n', num);

%get the number of lines of input file
[tmp,lines] = unix(sprintf('cut -f 1 %s | wc -w',gff_fname));
lines = int32(str2double(lines));

% for status of completeness
permille = int32(lines/1000);

% init struct genes
genes = init_genes(num);

% get all gene IDs
fid = fopen(gff_fname, 'r');
% write to log file
logfid = fopen('log.tmp','w');
% format for textscan of line
% 15 colums: <chromosome><origin><type><begin><end><score><strand><frame><'gene_id'><"><gID><'";'><'transcript_id'><"><tID><'";'>
format = '%s %*s %s %u %u %f %c %c %*s %*c %[^"] %*s %*s %*c %[^"] %*s';
% keeping only 9: <contig><type><begin><end><score><strand><frame><gID><tID>

% chromosome number
chr_num=-1;
% gene position
gpos=1;
% counting number of read lines
count=1;
% counting number of missing items
noutrb=0;
noutra=0;
nogid=0;
notid=0;
% instantiating a suffix tree for fast lookup of already contained geneIDs
% [0:9] is the alphabet which all the IDs are from (for they are numbers)
% sign '_' added to alphabet to distinguish gene duplications e.g. 0815, 0815_1, 0815_2, ...
gid_st=suffix_tree(['0':'9','_']);

while ~feof(fid)
  % parse line by line
  vec = textscan(fid, format, 1);
  
  % filter for contig
  if ~isempty(contig) && ~ismember(vec{1},contig)
    %fprintf(fid,'');
    count=count+1;
    continue;
  end

  % filter for sources
  if ~isempty(source) && ~ismember(vec{2}, source)
    count=count+1;
    continue;
  end

  % filter for types
  [relevant,tmp] = ismember(vec{2},SOFA_short_names);
  if ~relevant
    count=count+1;
    continue;
  end
  
  % counting number of missing gene IDs and skipping these lines
  if isempty(vec{8})
    nogid=nogid+1;
    count=count+1;
    continue;
  end

  % counting the number of missing transcript IDs and skipping these lines
  if isempty(vec{9})
    notid=notid+1;
    count=count+1;
    continue
  end
  
  switch char(vec{2}) ,
   case {'CDS','coding_exon'}
    etype = 3;
   case 'three_prime_UTR'
    etype = 4;
   case 'five_prime_UTR'
    etype = 2; 
   case 'exon'
    etype = 0;
   otherwise
    % skip line
    count=count+1;
    continue
  end

  % first adapt indices to stop-start=length scheme
  start = vec{3};
  stop  = vec{4};
  if strcmp(vec{6},'-')
    start=start-1;
  else
    stop=stop+1;
  end

  % determine strand
  if strcmp(vec{6},'+')
    plusstr=true;
  else
    plusstr=false;
  end

% insert geneID if not yet present

  gid=vec{8}{1:1}(3:end);

  if isempty(genes(gpos).name)
    % only for first insertion
    gid_st=gid_st.insert(vec{8}{1:1}(3:end),gpos); % insert in suffix tree
    genes(gpos).name=char(vec{8});
    genes(gpos).id=gpos;
    genes(gpos).chr=char(vec{1});
    chr_num=find(strcmp(char(vec{1}),genome_info.contig_names));
    genes(gpos).chr_num=chr_num;
    genes(gpos).strand=char(vec{6}); % now transcriptspecific
  elseif ~strcmp(genes(gpos).name,vec{8})
    % gene complete
    % will be set to true if something critical is missing
    ignore_gene=false; 
    if strcmp(genes(gpos).strand,'+')
      startfield='utr5_exons';
      stopfield='utr3_exons';
    else
      startfield='utr3_exons';
      stopfield='utr5_exons';
    end
    i=1;
    %if gpos==657
    %  warning('mal schaun was passiert')
    %end
    while i<=length(genes(gpos).transcripts)
      del_trans=false;
      if isempty(genes(gpos).(startfield){i})
        noutrb=noutrb+1;
        fprintf(logfid,'line %u: no UTR before gene %s\n', count, genes(gpos).name);
        del_trans=true;
      end
      if isempty(genes(gpos).(stopfield){i})
        noutra=noutra+1;
        fprintf(logfid,'line %u: no UTR after gene %s\n', count, genes(gpos).name);
        del_trans=true;
      end
      if del_trans
        % transcript will be deleted
        genes(gpos).transcripts(i)=[];
        genes(gpos).exons(i)=[];
        genes(gpos).cds_exons(i)=[];
        genes(gpos).utr5_exons(i)=[];
        genes(gpos).utr3_exons(i)=[];
      else
        i=i+1;
      end
    end
    % set borders
    for i=1:length(genes(gpos).transcripts)
      genes(gpos).start(end+1)=genes(gpos).exons{i}(1,1);
      genes(gpos).stop(end+1)=genes(gpos).exons{i}(end,2);
    end

    if isempty(genes(gpos).cds_exons)
      ignore_gene=true;
    else

      for i=1:length(genes(gpos).transcripts)
        if strcmp(genes(gpos).strand,'+')
          genes(gpos).tis(end+1)=genes(gpos).cds_exons{i}(1,1);
          genes(gpos).cdsStop(end+1)=genes(gpos).cds_exons{i}(end,2)-3;
          genes(gpos).tss(end+1)=genes(gpos).start(i);
          genes(gpos).cleave(end+1)=genes(gpos).stop(i);
        else
          genes(gpos).tis(end+1)=genes(gpos).cds_exons{i}(end,2);
          genes(gpos).cdsStop(end+1)=genes(gpos).cds_exons{i}(1,1)+3;
          genes(gpos).tss(end+1)=genes(gpos).stop(i);
          genes(gpos).cleave(end+1)=genes(gpos).start(i);
        end
      end
    end

    % check if suffix_tree gid_st already contains completed gene
    [tf,pos]=gid_st.contains(genes(gpos).name(3:end));
    %for k=1:length(pos)
    %  % for all equal gene_ids
    %  if strcmp(genes(pos(k)).strand,genes(gpos).strand)
    %    % on the same strand
    %    if genes(pos(k)).start > genes(gpos).stop || genes(gpos).start > genes(pos(k)).stop
    %      % the dont overlap => change name of gene
    %      % implicit by ignore_gene==false & tf
    %    else
    %      ignore_gene=true;
    %    end 
    %    fprintf(logfid,'line %u: gene %s already contained for the other strand. Skipping\n', count, gid);
    %  else
    %    % on different strands => change name of gene
    %    % implicit by ignore_gene==false & tf
    %  end
    %end
    if tf && ~ignore_gene
      % change name
      % determine number of prefices in gid_st=#(geneID duplicates)
      duplicates=gid_st.contains_prfx(genes(gpos).name(3:end));
      genes(gpos).name=strcat(genes(gpos).name,'_',int2str(duplicates)); % is correct for first ID has no number
    end
    % insert completed gene into suffix tree
    gid_st=gid_st.insert(genes(gpos).name(3:end),gpos);

    % check if chromosom has changed
    if ~strmatch(char(vec{1}), genes(gpos).chr)
      % uopdate chr_num;
      chr_num=find(strcmp(char(vec{1}),genome_info.contig_names));
    end

    % initialize next gene
    if ignore_gene
      % gene will be ignored so delete entries
      fprintf(logfid,'line %u: Gene %s will be overwritten', count, gid);
      genes(gpos).transcripts = [];
      genes(gpos).exons = [];
      genes(gpos).cds_exons = [];
      genes(gpos).utr5_exons = [];
      genes(gpos).utr3_exons = [];
    else
      gpos=gpos+1;      
    end
    genes(gpos).name=char(vec{8});
    genes(gpos).id=gpos;
    genes(gpos).chr=char(vec{1});
    genes(gpos).chr_num=chr_num;
    genes(gpos).strand=char(vec{6});
  end
 
  % search transcript

  % first fill
  if isempty(genes(gpos).transcripts)
    genes(gpos).transcripts = cellstr(vec{9});
    genes(gpos).exons{end+1}=zeros(0,3);
    genes(gpos).cds_exons{end+1}=zeros(0,3);
    genes(gpos).utr5_exons{end+1}=zeros(0,3);
    genes(gpos).utr3_exons{end+1}=zeros(0,3);
    % transcriptspecific strand
    %genes(gpos).strand{end+1}=char(vec{6});
    tpos=1;
  else
    % try to find transcript
    [found,tpos]=ismember(vec{9}, genes(gpos).transcripts);
    % not yet contained
    if ~found
      genes(gpos).transcripts(end+1)=cellstr(vec{9});
      tpos=length(genes(gpos).transcripts);
      genes(gpos).exons{end+1}=zeros(0,3);
      genes(gpos).cds_exons{end+1}=zeros(0,3);
      genes(gpos).utr5_exons{end+1}=zeros(0,3);
      genes(gpos).utr3_exons{end+1}=zeros(0,3);
      % transcriptspecific strand
      %genes(gpos).strand{end+1}=char(vec{6});
    end
  end
  % add exons

  % adding to exons
  if etype==0
    genes(gpos).exons{tpos}(end+1,:)=[start, stop, etype];
  end

  % adding depending on type
  % determine field name
  switch etype
   case 2 % 5'UTR
    field='utr5_exons';
   case 4 % 3'UTR
    field='utr3_exons';
   case 3 % CDS
    field='cds_exons';
   case 0 % exon (assumption cds always before exon)
    % now we have to look if it is already in cds_exons
    if isempty(genes(gpos).cds_exons{tpos})
      sfound=false;
      efound=false;
    else
      sfound=(start==genes(gpos).cds_exons{tpos}(end,1));
      efound=(stop==genes(gpos).cds_exons{tpos}(end,2));
    end
    if sfound && efound
      % its a cds_exon and already in
      count=count+1;
      continue;
    elseif ~sfound && ~efound
      % its not a cds_exon so add it as a utr one
      if isempty(genes(gpos).cds_exons{tpos})
        % UTR before CDS
        if plusstr
          field='utr5_exons';
          etype=2;
        else
          field='utr3_exons';
          etype=4;
        end
      else
        % UTR after CDS
        if plusstr
          field='utr3_exons';
          etype=4;
        else
          field='utr5_exons';
          etype=2;
        end
      end

    % cases in which exon has to be divided
    % 'field' is no use, continue required
    elseif efound

      % its first part is UTR second CDS
      if plusstr
        genes(gpos).utr5_exons{tpos}(end+1,:)=[start, genes(gpos).cds_exons{tpos}(end,1), 2];
      else
        genes(gpos).utr3_exons{tpos}(end+1,:)=[start, genes(gpos).cds_exons{tpos}(end,1), 2];
      end
      count=count+1;
      continue;

    elseif sfound

      % its first part is CDS second UTR
      if plusstr
        genes(gpos).utr3_exons{tpos}(end+1,:)=[genes(gpos).cds_exons{tpos}(end,2), stop, 2];
      else
        genes(gpos).utr5_exons{tpos}(end+1,:)=[genes(gpos).cds_exons{tpos}(end,2), stop, 2];
      end
      count=count+1;
      continue;

    else
      warning('no case matches')
    end
   otherwise
    warning('inappropriate structure of infile data');
  end

  % Adding exon depending on field
  if plusstr
    genes(gpos).(field){tpos}(end+1,:)=[start, stop, etype];
  else
    genes(gpos).(field){tpos}(end+1,:)=[start, stop, etype];
  end
  if mod(count,permille)==0
    p_value= double(count/permille)/10;
    fprintf('\r%3.1f %% of infile parsed', p_value);
  end
  %% for testing
  %if mod(gpos,1001)==0
  %  gpos=gpos-1;
  %  fprintf('\nending testing\n');
  %  fprintf('count=%u\n',count);
  %  fprintf('gpos=%u\n',gpos);
  %  %for last gene
  %  genes(gpos+1:end)=[];
  %  if ~isempty(save_fname)
  %    save(save_fname, 'genes');
  %  end
  %  return
  %end
  count=count+1;
end
fclose(fid);
% set for the last gene
% will be set to true if something critical is missing
ignore_gene=false; 

if strcmp(genes(gpos).strand,'+')
  startfield='utr5_exons';
  stopfield='utr3_exons';
else
  startfield='utr3_exons';
  stopfield='utr5_exons';
end
i=1;
while i<=length(genes(gpos).transcripts)
  del_trans=false;
  if isempty(genes(gpos).(startfield){i})
    noutrb=noutrb+1;
    fprintf(logfid,'line %u: no UTR before gene %s\n', count, genes(gpos).name);
    del_trans=true;
  end
  if isempty(genes(gpos).(stopfield){i})
    noutra=noutra+1;
    fprintf(logfid,'line %u: no UTR after gene %s\n', count, genes(gpos).name);
    del_trans=true;
  end
  if del_trans
    % transcript will be deleted
    genes(gpos).transcripts(i)=[];
    genes(gpos).exons(i)=[];
    genes(gpos).cds_exons(i)=[];
    genes(gpos).utr5_exons(i)=[];
    genes(gpos).utr3_exons(i)=[];
  else
    i=i+1;
  end
end

% set borders
for i=1:length(genes(gpos).transcripts)
  genes(gpos).start(end+1)=genes(gpos).exons{i}(1,1);
  genes(gpos).stop(end+1)=genes(gpos).exons{i}(end,2);
end

if isempty(genes(gpos).cds_exons)
  ignore_gene=true;
else

  for i=1:length(genes(gpos).transcripts)
    if strcmp(genes(gpos).strand,'+')
      genes(gpos).tis(end+1)=genes(gpos).cds_exons{i}(1,1);
      genes(gpos).cdsStop(end+1)=genes(gpos).cds_exons{i}(end,2)-3;
      genes(gpos).tss(end+1)=genes(gpos).start(i);
      genes(gpos).cleave(end+1)=genes(gpos).stop(i);
    else
      genes(gpos).tis(end+1)=genes(gpos).cds_exons{i}(end,2);
      genes(gpos).cdsStop(end+1)=genes(gpos).cds_exons{i}(1,1)+3;
      genes(gpos).tss(end+1)=genes(gpos).stop(i);
      genes(gpos).cleave(end+1)=genes(gpos).start(i);
    end
  end
end

% check if suffix_tree gid_st already contains completed gene
[tf,pos]=gid_st.contains(genes(gpos).name(3:end));
for k=1:length(pos)
  % for all equal gene_ids
  if strcmp(genes(pos(k)).strand,genes(gpos).strand)
    % on the same strand
    if genes(pos(k)).start > genes(gpos).stop || genes(gpos).start > genes(pos(k)).stop
      % the dont overlap => change name of gene
      % implicit by ignore_gene==false & tf
    else
      ignore_gene=true;
    end 
    fprintf(logfid,'line %u: gene %s already contained for the other strand. Skipping\n', count, gid);
  else
    % on different strands => change name of gene
    % implicit by ignore_gene==false & tf
  end
end
if tf && ~ignore_gene
  % change name
  % determine number of prefices in gid_st=#(geneID duplicates)
  duplicates=gid_st.contains_prfx(genes(gpos).name(3:end));
  genes(gpos).name=strcat(genes(gpos).name,'_',int2str(duplicates)); % is correct for first ID has no number
end
% insert completed gene into suffix tree
gid_st=gid_st.insert(genes(gpos).name(3:end),gpos);
if ignore_gene
  % gene will be ignored so delete entries
  fprintf(logfid,'line %u: Gene %s will be overwritten', count, gid);
  genes(gpos).transcripts = [];
  genes(gpos).exons = [];
  genes(gpos).cds_exons = [];
  genes(gpos).utr5_exons = [];
  genes(gpos).utr3_exons = [];
else
  gpos=gpos+1;
end
fprintf('\nnumber of genes in infile: %u', num);
fprintf('\nnumber of parsed genes: %u', gpos);
fprintf('\nnumber of accepted lines: %u', count);
fprintf('\nno UTR preceding %u transcripts', noutrb);
fprintf('\nno UTR succeeding %u genes', noutra);
fprintf('\n%u genes without ID', nogid);
fprintf('\n%u transcripts without ID\n', notid);
% cut away unused entries
genes(gpos:end)=[];
if ~isempty(save_fname)
   save(save_fname, 'genes');
end
fclose(logfid);
end
