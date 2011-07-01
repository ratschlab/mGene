function genometool(fasta_fname, fasta_dir, information_file, base_dir) ;
% genometool(fasta_fname, fasta_dir, information_file, genome_config_dir) ;

[ret, timedate] = unix('date') ;
assert(ret==0) ;
timedate(timedate==sprintf('\n')) = [] ;
fprintf('started genometool at %s\n', timedate) ;

fprintf('Results are written to: %s\n', information_file) ;

% find internal files
[fasta_fname, fasta_dir] = find_internal_files(fasta_fname, fasta_dir) ;

mkdir(base_dir) ;

if exist('information_file')
  fd = fopen(information_file, 'w+') ;
  if fd<1, error('could not open file %s', information_file); end
  fprintf(fd,'---------------------------------------------------- \n');
  fprintf(fd,'GenomeTool started %s\n', timedate) ;
  fprintf(fd,'---------------------------------------------------- \n\n');
  fclose(fd)   
end

genome_config = [base_dir '/genome.config'] ;

fprintf('Using fasta file from %s \n\n',fasta_fname) ;

if ~isempty(fasta_dir) && fexist(sprintf('%s/genome.config', fasta_dir)),
  fprintf('Creating a link to existing genome information object\n') ;

  genome_info = init_genome(sprintf('%s/genome.config', fasta_dir)) ;
  ret = unix(sprintf('ln -s %s %s; ln -s %s %s', sprintf('%s/genome.config', fasta_dir), genome_config, [fasta_dir '/genome'], [base_dir '/genome'])) ;
  assert(ret==0) ;

  if exist('information_file') %~isempty(fd),
    fd = fopen(information_file, 'a+') ;

    fprintf(fd, 'Reading FASTA File : %s\n', fasta_fname) ; 
    fprintf(fd, 'Genome properties:\n') ;
    fprintf(fd, ' * %i contigs\n', length(genome_info.contig_names)) ;
  
    len = 0 ;
    for i=1:length(genome_info.flat_fnames),
      if ~(fexist(genome_info.flat_fnames{i})==1), 
	error('could not find file %s', genome_info.flat_fnames{i});
      end
      d=dir(genome_info.flat_fnames{i}) ;
      assert(length(d)==1) ;
      len = len + d.bytes ;
    end 
    fprintf(fd, ' * %ikb total length\n\n', round(len/1024)) ;
    fprintf(fd, 'Contig list:\n') ;
    for i=1:length(genome_info.contig_names),
      fprintf(fd, '  %s\n', genome_info.contig_names{i}) ;
    end 
    fclose(fd) ;
  end
else
  python_script=sprintf('%s/tools/make_gio.py',get_base_dir());
  % python_script=sprintf('make_gio.py');
  python_call = sprintf('python %s %s %s %s',python_script, fasta_fname, base_dir, information_file)
  ret=unix(python_call);
  assert(ret==0) ;
end ;

if exist('information_file'),
  fd = fopen(information_file, 'a+') ;
  if fd<1, error('could not open file %s', information_file); end
  [ret, timedate] = unix('date') ;
  assert(ret==0) ;
  timedate(timedate==sprintf('\n')) = [] ;
  fprintf(fd,'\n---------------------------------------------------- \n');
  fprintf(fd,'finished GenomeTool at %s\n', timedate) ;
  fprintf(fd,'---------------------------------------------------- \n');
  fclose(fd) ;

  ret=unix(sprintf('~/bin/cleanup_output.sh %s', information_file)) ;
  assert(ret==0) ;
end ;
