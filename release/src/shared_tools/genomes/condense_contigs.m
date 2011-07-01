function condense_contigs(old_config_fname,target_dir,Nspacer_length,target_contig_num)

% condense_contigs(fname,target_dir,Nspacer_length,target_contig_num)
% 
% this function should be used to concatinate contigs. It generates new
% fasta and fa files, as well as a new config file in the target_dir. In the fasta and fa
% files original contigs are separated by a stretch of Ns of length
% Nspacer_length (default: 25000; for genefinding this should be longer
% than the max. intron length, to avoid gene predictions across contigs)
% The number of resulting merged contigs can be specified by
% target_contig_num (default: 50). The original contigs are sorted such
% that the resulting contigs are of approx. equal lengths.
% The mapping between original and target contigs is specified by the
% file condense_map.txt in the target directory.
%  
% 
% INPUT  - old_config_fname
%        - target_dir 
%        - Nspacer_length (default: 25000)
%        - target_contig_num (default: 50)

if nargin<3
  Nspacer_length = 25000 ;
end
if nargin<4
  target_contig_num = 50 ;
end

rand('seed',17562354);

target_fname = sprintf('%s/genome.config',target_dir);
Ns_fname = sprintf('%s/N_map.gff3',target_dir);
% target fnames and directories 
target_map_fname = [target_dir '/condense_map.txt'] ;
target_genome_dir = [target_dir '/genome'] ;
mkdir(target_dir) ;
mkdir(target_genome_dir) ;

% collect information
genome_info = init_genome(old_config_fname) ;

contig_sizes = zeros(1,length(genome_info.flat_fnames)) ;
for i=1:length(genome_info.flat_fnames)
  d=dir(genome_info.flat_fnames{i}) ;
  contig_sizes(i) = d.bytes ;
end ;
[tmp,sidx]=sort(contig_sizes, 'descend') ;

target_contig_sizes=zeros(1,target_contig_num) ;
target_contig_sets={} ;
target_contig_sets{target_contig_num}=[] ;

% sort into new contigs
for i=1:length(genome_info.flat_fnames),
  if contig_sizes(sidx(i))==0, continue ; end ;
  [tmp,idx]=min(target_contig_sizes) ;
  target_contig_sizes(idx) = target_contig_sizes(idx) + contig_sizes(sidx(i)) ;
  target_contig_sets{idx}(end+1) = sidx(i) ;
  if length(target_contig_sets{idx})>1,
    target_contig_sizes(idx) = target_contig_sizes(idx) + Nspacer_length ;
  end ;
end 

rmidx=find(target_contig_sizes==0) ;
target_contig_sizes(rmidx) = [] ;
target_contig_sets(rmidx) = [] ;

% create new genome_info
genome_info_new = genome_info ;
genome_info_new.contig_names={} ;
genome_info_new.fasta_fnames={} ;
genome_info_new.flat_fnames={} ;
genome_info_new.basedir = target_dir ;

% write files
fd=fopen(target_map_fname, 'w+') ;
for i=1:length(target_contig_sets),
  idx=randperm(length(target_contig_sets{i})) ;
  target_contig_sets{i} = target_contig_sets{i}(idx) ;

  contig_name = sprintf('contig%i', i-1) ;
  genome_info_new.contig_names{end+1} = contig_name ;
  genome_info_new.fasta_fnames{end+1} = [target_genome_dir sprintf('/%s.fa', contig_name)] ;
  genome_info_new.flat_fnames{end+1} = [target_genome_dir sprintf('/%s.flat', contig_name)] ;
  
  target_flat_fd = fopen(genome_info_new.flat_fnames{end}, 'w+') ;
  target_flat_pos = 0 ;
  sets=target_contig_sets{i} ;
  for j=1:length(sets),
    source_flat_fd = fopen(genome_info.flat_fnames{sets(j)}, 'r') ;
    seq=char(fread(source_flat_fd,inf)) ;
    fclose(source_flat_fd) ;
    fwrite(target_flat_fd, seq) ;
    fprintf(fd, '%s\t%i\t%i\t%s\t%i\t%i\n', contig_name, target_flat_pos+1, target_flat_pos+length(seq), genome_info.contig_names{sets(j)}, 1, length(seq)) ;
    target_flat_pos = target_flat_pos+length(seq) ;
    if j~=length(sets),
      seq = 'N'*ones(1,Nspacer_length) ;
      fwrite(target_flat_fd, seq) ;
      fprintf(fd, '%s\t%i\t%i\tNSPACER\t%i\t%i\n', contig_name, target_flat_pos+1, target_flat_pos+length(seq), 1, Nspacer_length) ;
      target_flat_pos = target_flat_pos+Nspacer_length ;
    end ;
  end ;
  fclose(target_flat_fd) ;
  target_flat_fd = fopen(genome_info_new.flat_fnames{end}, 'r') ;
  seq=char(fread(target_flat_fd,inf)) ;
  fclose(target_flat_fd) ;

  target_fasta_fd = fopen(genome_info_new.fasta_fnames{end}, 'w+') ;
  write_fasta(target_fasta_fd, contig_name, seq) ;
  fclose(target_fasta_fd) ;
end 
fclose(fd) ;

% replace the contig lines in the genome config file
fd=fopen(old_config_fname, 'r') ;
target_fd = fopen(target_fname, 'w+') ;
while ~feof(fd),
  line = fgetl(fd) ;
  if ~ischar(line), break ; end ;
  % if isempty(line), continue ; end ;
  if length(line)>8 & isequal(line(1:7),'BASEDIR')
    fprintf(target_fd, 'BASEDIR %s\n', target_dir) ;
  elseif length(line)>7 & isequal(line(1:7),'CONTIGS'),
    for i=1:str2num(line(8:end))
      l2=fgetl(fd) ; 
    end 
    fprintf(target_fd, 'CONTIGS %i\n',length(genome_info_new.flat_fnames))
    for i=1:length(genome_info_new.flat_fnames),
      fprintf(target_fd, '%s\t%s\t%s\n', genome_info_new.contig_names{i}, ...
              strrep(strrep(genome_info_new.flat_fnames{i}, '//', '/'), target_dir, ''), genome_info_new.fasta_fnames{i}) ;
    end ;
  else
    fprintf(target_fd, '%s\n', line) ;
  end ;
end ;
fclose(fd);
fclose(target_fd);

generate_listofNs(target_map_fname,Ns_fname) ;

link_name = sprintf('%s/genome_uncondensed.config',target_dir);
unix(sprintf('ln -s %s %s', old_config_fname,link_name));

