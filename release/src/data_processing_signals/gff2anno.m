function gff2anno(genome_dir, gff_fname, annotation_fname, annotation_dir, version, seqid, source, level, ignore_missing_referenced_entries, run_locally, converter) ;
%        gff2anno(genome_dir, gff_fname, annotation_fname, annotation_dir, version, seqid, source, level, ignore_missing_referenced_entries, run_locally, converter) ;


[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
nice_mkdir(annotation_dir)
fd=fopen(annotation_fname, 'w+') ;
if fd<1, error('could not open file %s', annotation_fname); end

fprintf(fd,'-------------------------------------------------- \n') ;
fprintf(fd,'GFF2Anno started %s\n', timedate) ;
fprintf(fd,'-------------------------------------------------- \n\n') ;
fprintf('-------------------------------------------------- \n') ;
fprintf('GFF2Anno started %s\n', timedate) ;
fprintf('-------------------------------------------------- \n\n') ;

if ~exist('converter', 'var'),
  converter='';
end ;

if ~exist('run_locally')
  run_locally =-1;
end

if isequal(source, '-'),
  source='' ;
end ;
if isequal(seqid, '-'),
  seqid='' ;
end;
if isequal(level, '-'),
  level='' ;
end;
if ~exist('ignore_missing_referenced_entries'),
  ignore_missing_referenced_entries = 1 ;
end ;


fprintf(fd, 'Reading GFF File : %s\n\n',gff_fname) ; 
fprintf('Reading GFF File : %s\n\n',gff_fname) ; 

do_correct_tis_stop = 0 ;
if isempty(converter) && isequal(version, 'ensembl-gtf')
  converter = 'gtf-to-gff3' ;
  do_correct_tis_stop = 1 ; % ensembl GTF files do not often show
                            % the consensus sites for TIS and STOP
end ;

if ~isempty(converter),
  if any(converter == ',')
    converter = separate(converter, ','); 
  else
    converter = {converter} ;
  end ;

  if ismember('gtf-to-gff3', converter)
    fprintf(fd, 'Running GTF-to-GFF3 converter before parsing the file\n') ;
    fprintf('Running GTF-to-GFF3 converter before parsing the file\n') ;

    tmpname = sprintf('%s/gtf.gff3', annotation_dir) ;
    fprintf('Writing GFF3 file to %s\n', tmpname) ;

    [engine,tmp,tmp,mccdir]=  determine_engine;
    home=getenv('HOME');
    home = deblank(home);
    cmp_fct = sprintf('%s/mgene_galaxy/gtf2gff/run_gtf2gff_helper.sh',home) ;

    if fexist(cmp_fct) && ~isequal(engine, 'matlab') ;
      fprintf('using precompiled code (system)\n')
      [ret] = system(sprintf('%s %s %s %s',...
			     cmp_fct, mccdir, gff_fname, tmpname)) ;
       assert(ret==0) ;
    else
      fprintf('not using precompiled code\n')
      gtf2gff(gff_fname, tmpname) ;
    end

    gff_fname = tmpname ;
  end ;

  if ismember('correct-tis-stop', converter),
    do_correct_tis_stop = 1 ;
  end ;
end ;

disp('-------------------------------------') ;
disp('Step 1: Loading genome information...') ;
disp('-------------------------------------') ;

genome_config = [genome_dir '/genome.config'] ;
genome_info = init_genome(genome_config);



disp('------------------------------------') ;
disp('Step 2: Setting up datastructures...') ;
disp('------------------------------------') ;

fprintf('Input parameters:\n') ;
fprintf('version: "%s"\n', version) ;
fprintf('source : "%s"\n', source) ;
fprintf('seqid  : "%s"\n', seqid) ;
fprintf('level  : "%s"\n', level) ;

if isequal(version, 'generic')
  seqid_names = {};
  source_names = {} ;
  level1_names = {'gene'} ;
  level2_names = {'mRNA'} ;
  level3_names = {'CDS','five_prime_UTR','three_prime_UTR'} ;
  level3_coding_names = {'CDS'} ;
elseif isequal(version, 'wormbase')
  seqid_names = {};
  source_names = {'Coding_transcript'} ;
  level1_names = {'gene'} ;
  level2_names = {'mRNA'} ;
  level3_names = {'CDS','five_prime_UTR','three_prime_UTR'} ;
  level3_coding_names = {'CDS'} ;
elseif isequal(version, 'flybase')
  seqid_names = {};
  source_names = {'FlyBase','flybase'} ;
  level1_names = {'gene'} ;
  level2_names = {'mRNA'} ;
  level3_names = {'CDS','five_prime_UTR','three_prime_UTR'} ;
  level3_coding_names = {'CDS'} ;
elseif isequal(version, 'tair')
  seqid_names = {};
  source_names = {'TAIR8'} ;
  level1_names = {'gene'} ;
  level2_names = {'mRNA'} ;
  level3_names = {'CDS','five_prime_UTR','three_prime_UTR'} ;
  level3_coding_names = {'CDS'} ;
elseif isequal(version, 'ensembl-gtf')
  seqid_names = {};
  source_names = {'protein_coding'} ;
  level1_names = {'gene'} ;
  level2_names = {'mRNA'} ;
  level3_names = {'CDS','exon','five_prime_UTR','three_prime_UTR'} ;
  level3_coding_names = {'CDS'} ;
elseif isequal(version, 'other') 
  seqid_names = {};
  if ~isempty(seqid) && ismember(',',seqid)
    seqid_names = separate(seqid,','); 
  elseif ~isempty(seqid)
    seqid_names{1} = seqid;
  end ;
  
  source_names = {};
  if ~isempty(source) &&ismember(',', source)
    source_names = separate(source,','); 
  elseif ~isempty(source)
  source_names{1} = source;
  end
  
  levels = separate(level,':') ; 
  assert(isequal(levels{1},'level'))
  level1_names = {};
  if ~isempty(levels{2})&&ismember(',',levels{2})
    level1_names = separate(levels{2},','); 
  elseif ~isempty(levels{2})
    level1_names = levels(2);
  end
  
  level2_names = {};
  if ~isempty(levels{3})&&ismember(',',levels{3})
    level2_names = separate(levels{3},','); 
  elseif ~isempty(levels{3})
    level2_names = levels(3);
  end
  level3_names = {};
  level3_coding_names = {};
  if ~isempty(levels{4})&&ismember(',',levels{4})
    level3_names = separate(levels{4},','); 
    level3_coding_names = separate(levels{4},','); ;
  elseif ~isempty(levels{4})
    level3_names = levels(4);
    level3_coding_names = levels(4);
  end
  if isequal(levels{5}, 'None'),
    levels{5}='' ;
  end ;
  if ~isempty(levels{5}) && ismember(',',levels{5})
    names = separate(levels{5},','); 
    level3_names = {level3_names{:}, names{:}};
  elseif ~isempty(levels{5})
    level3_names = {level3_names{:}, levels{5}};
  end
else
  error('unknown version %s', version)
end

if ~isempty(seqid_names)
  fprintf(fd, 'Parsing following sequids: \n') ; 
  for i=1:length(seqid_names)
	  fprintf(fd, '\t %s\n',seqid_names{i}) ; 
  end
end ;

if ~isempty(source_names)
  fprintf(fd, 'Parsing following sources: \n') ; 
  for i=1:length(source_names)
  	fprintf(fd, '\t %s\n',source_names{i}) ; 
  end
end

fprintf(fd, 'Parsing following types (level 1): \n') ; 
for i=1:length(level1_names)
	fprintf(fd, '\t %s\n',level1_names{i}) ; 
end
fprintf(fd, 'Parsing following types (level 2): \n') ; 
for i=1:length(level2_names)
	fprintf(fd, '\t %s\n',level2_names{i}) ; 
end
fprintf(fd, 'Parsing following types (level 3): \n') ; 
for i=1:length(level3_names)
	fprintf(fd, '\t %s\n',level3_names{i}) ; 
end
fprintf(fd, 'Following types are considered coding exons: \n') ; 
for i=1:length(level3_coding_names)
	fprintf(fd, '\t %s\n\n',level3_coding_names{i}) ; 
end
fclose(fd) ;

TYPE.short_names = {level1_names{:} level2_names{:} level3_names{:}} ;
TYPE.IDs = TYPE.short_names ;
TYPE.hierachies = [ones(1,length(level1_names)) 2*ones(1,length(level2_names)) 3*ones(1,length(level3_names))] ;

fprintf('Done.\n\n') ;

nofgenes=0 ;
for i=1:length(level1_names),
	[ret, nofgenes_] = unix(sprintf('grep %s %s | wc -l', level1_names{i}, gff_fname)) ;
  assert(ret==0)
  nofgenes_ = str2num(nofgenes_);
  nofgenes = nofgenes + nofgenes_ ;
end ;

fprintf('Found %i level1 objects (genes, etc.) \n', nofgenes) ;

genes = init_genes(max(nofgenes,1));
[ret, seq_ids] = unix(sprintf('sed -e "1d" %s | grep -v "^#" | grep -v "^$" | cut -f 1 | sort -u', gff_fname)); 
assert(ret==0) ;

if isempty(seqid_names)
  contig_names = separate(seq_ids, sprintf('\n')) ;
  if ~isempty(contig_names) && isempty(contig_names{end})
    contig_names(end)=[] ;
  end ;
else
  contig_names = seqid_names;
end

% what is this useful for?
% in any case this is extremely slow for bigger files
if 0,
  N=0;
  tmp_fname = tempname;
  for contig=1:length(contig_names)
    cname = contig_names{contig};
    unix(sprintf('grep "^%s\t" %s > %s', cname,gff_fname,tmp_fname));
    [tmp,num_seqs] = unix(sprintf('sed -e "1d" %s | cut -f 1 | wc -l',tmp_fname));
    N = N+str2num( num_seqs);
  end
  if N==0
    error('GFF file does not contain specifed seqids')
  end
  unix(sprintf('rm -f %i', tmp_fname));
end ;

if run_locally ==-1
  run_locally = rproc_policy('gff2anno',genome_info,[nofgenes length(contig_names)]) ;
end

Data.source = source_names;
Data.gff_fname = gff_fname;
Data.TYPE = TYPE;
Data.genome_config = genome_config;
Data.level3_coding_names= level3_coding_names;
Data.do_correct_tis_stop = do_correct_tis_stop ;

%run_locally = 1 ; 

% split by contig
if ~run_locally
  fprintf('parsing by contig\n')
  for contig=1:length(contig_names)
    Data.cname = contig_names(contig);
    Data.save_fname = sprintf('%s/%s_', annotation_dir, Data.cname{:});
    Data.log_fname  = sprintf('%s/%s_log', annotation_dir, Data.cname{:});
	options.start_dir = get_base_dir();
	options.addpaths =  {fileparts(which('start_gene_predict_rgasp'))};

    [mem_req, time_req, opts] = rproc_memtime_policy('gff2anno_helper',  0, options) ;

    jobinfo(contig)=rproc('gff2anno_helper', Data, mem_req, opts, time_req);
  end
  jobinfo = rproc_wait(jobinfo, 30, 1, -1);
  len = 0 ;
  for contig=1:length(contig_names)
    save_fname = sprintf('%s/%s_',annotation_dir,contig_names{contig});
    l = load([save_fname 'genes.mat'],'genes'); 
    genes(len+1:len+length(l.genes))=l.genes;
    len = len + length(l.genes);
    l = load([save_fname 'FAULTS1.mat'],'FAULTS1'); 
    FAULTS1(contig) = l.FAULTS1;
    l = load([save_fname 'FAULTS2.mat'],'FAULTS2'); 
    FAULTS2(contig) = l.FAULTS2;
    clear l
    log_fname  = sprintf('%s/%s_log', annotation_dir, Data.cname{:});
    unix(sprintf('cat %s >> %s', log_fname, annotation_fname)) ;
    unix(sprintf('rm %sgenes.mat',save_fname));
    unix(sprintf('rm %sFAULTS1.mat',save_fname));
    unix(sprintf('rm %sFAULTS2.mat',save_fname));
    genes = genes(1:len) ;
  end
else
  Data.cname = contig_names;
  Data.save_fname = [annotation_dir '/all_'];
  Data.log_fname  = [annotation_dir '/all_log'];

  gff2anno_helper(Data) ;

  load([Data.save_fname 'genes.mat'],'genes'); 
  load([Data.save_fname 'FAULTS1.mat'],'FAULTS1'); 
  load([Data.save_fname 'FAULTS2.mat'],'FAULTS2'); 

  unix(sprintf('cat %s >> %s', Data.log_fname, annotation_fname )) ;
  unix(sprintf('rm %s', Data.log_fname));
end
for i=1:length(genes)
  genes(i).id = i;
end

fprintf('Done.\n\n') ;
%annotation_fname


fd=fopen(annotation_fname, 'a') ;
if fd<1, error('could not open file %s', annotation_fname); end
fprintf(fd, '\nAnnotation with %i genes ', length(genes)) ;
% add some more information here

disp('------------------------------------') ;
disp('Step 7: Saving genome annotation ...') ;
disp('------------------------------------') ;

% unix(sprintf('mkdir -p %s', annotation_dir)) ;
save_struct(sprintf('%s/genes.mat', annotation_dir), genes, 'genes') ;
save_struct(sprintf('%s/FAULTS1.mat', annotation_dir),FAULTS1, 'FAULTS1') ;
save_struct(sprintf('%s/FAULTS2.mat', annotation_dir),FAULTS2, 'FAULTS2') ;

% save(sprintf('%s/info.mat', annotation_dir), 'genome_info', 'genome_config') ;
fprintf(fd, '(saved to %s/genes.mat)\n', annotation_dir) ;

[ret, timedate] = unix('date') ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf(fd,'\n-------------------------------------------------- \n') ;
fprintf(fd,'GFF2Anno finished at %s\n', timedate) ;
fprintf(fd,'-------------------------------------------------- \n') ;
fprintf('\n-------------------------------------------------- \n') ;
fprintf('GFF2Anno finished at %s\n', timedate) ; 
fprintf('--------------------------------------------------- \n') ;

fclose(fd) ;

ret=unix(sprintf('~/bin/cleanup_output.sh %s', annotation_fname)) ;
assert(ret==0) ;

disp('Done.') ;


