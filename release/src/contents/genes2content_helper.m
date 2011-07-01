function genes2content_helper(genes, genome_config, cpf_fname, cpf_dir, content_name, vs_content_names,intergenic_dist_to_genes) ;
% genes2content_heper(genes, genome_config, cpf_fname, cpf_dir, content_name,  vs_content_names, intergenic_dist_to_genes) ;

disp('-------------------------------------------------------') ;
disp('Step 5: Extracting labeled segments from the genome...') ;
disp('-------------------------------------------------------') ;

if ~isempty(cpf_dir),
  galaxy_PAR.dir = cpf_dir ;
else
  galaxy_PAR.dir = tempname ;
end ;
galaxy_PAR.genome_config = genome_config ;
PAR = create_PAR('exp_galaxy', galaxy_PAR);
Content = PAR.Contents.(content_name);
segment_ids = PAR.model.segments;
ids = [];
for i=fieldnames(segment_ids)'
  ids = [ids segment_ids.(i{1})]; 
end
% for the moment contents must be one of segment types. 
if isfield( segment_ids,content_name)
  which_segments{1} = segment_ids.(content_name);
elseif isequal(content_name,'frame0')
  which_segments{1} = segment_ids.cds_exon;
end
if isempty(vs_content_names)
    which_segments{2} = setdiff(ids,which_segments{1});% one versus rest
else
  which_segments{2} = segment_ids.(vs_content_names);
end

if ~isempty(cpf_dir),
  out_dir = [cpf_dir '/tmp'] ;
else
  out_dir = tempname ;
end ;
unix(sprintf('mkdir -p %s', out_dir)) ;

if isempty(intergenic_dist_to_genes)
  intergenic_dist_to_genes = [100 20000];
end

if isempty(strfind(content_name,'frame0'))
  use_frame=0;
else
  use_frame=1;
end

genome_info = init_genome(genome_config) ;
contig_names = genome_info.contig_names ;
for c=1:length(contig_names)
  for s='+-',
    cs_idx = find(ismember({genes.chr}, contig_names{c}) & [genes.strand]==s) ;
    %if isempty(cs_idx)
    %  continue
    %end
    cs_genes = genes(cs_idx) ;
    P_ = struct ;
    P_.genes = cs_genes;
    P_.segment_ids = segment_ids;
    cs_genes = gen_gene_segmentation(P_);
    clear P_
    
    if isempty(cs_genes),
      segments=[] ;
      label=[] ;
    else
      % generate positive and negative candidate segments;
      P_ = struct ;
      P_.genes = cs_genes;
      P_.intergenic_dist_to_genes = intergenic_dist_to_genes;
      P_.use_frame = use_frame;
      P_.which_segments = which_segments;
      P_.segment_ids = segment_ids;
      P_.genome_info = genome_info;
	
	  % Content label function is 
	  % usually the function get_cand_contents
      [segments, label] = feval(Content.label_fct, P_);
      clear P_ ;
    end ;

    % save them to a file
    filename_prefix = sprintf('%s/contig_%i%c', out_dir, c, s) ;
    fn_in=[filename_prefix '_label.cpf'];
    % unix(sprintf('rm -f %s', fn_in)) ;
    save_predictions_as_contpred(filename_prefix, segments, label, 'label', content_name, contig_names{c}, s, 0) ;
    save_predictions_as_contpred_bin(filename_prefix, segments, label, 'label', content_name, contig_names{c}, s) ;
  end ;
end ;
fprintf('Done.\n\n') ;

disp('----------------------------------') ;
disp('Step 6: Generating output files...') ;
disp('----------------------------------') ;
unix(sprintf('rm -f %s', cpf_fname)) ;
cnt = 0 ;
for c = 1:length(contig_names)
  fprintf('Processing outputs for contig %s\n', contig_names{c}) ;
  for s = '+-'
    P.chrom = c;
    P.strand = s;
    filename = sprintf('%s/contig_%i%s', out_dir, c, s);
    fn_in = [filename '_label.cpf'];
    ret = unix(['cat ' fn_in ' >> ' cpf_fname ' && rm -f ' fn_in]) ;
    assert(ret==0) ;

    fn_mat = sprintf('%s/contig_%i%s_label_cpf.mat', out_dir, c, s);
    fnall_mat = sprintf('%s/label_cpf.mat', cpf_dir);
    if ~fexist(fn_mat)
      fprintf('binary content prediction file %s not found\n', fn_mat) ;
      continue
    end
    CPF=load(fn_mat);
    fieldnames(CPF);

    save_append(fnall_mat, 1, sprintf('CPF_%i', cnt), CPF) ;
    cnt = cnt + 1 ;
    ret = unix(['rm -f ' fn_mat]) ;
    assert(ret==0) ;
  end ;
end ;
if length(cpf_fname)>=3 && isequal(cpf_fname(end-2:end), '.gz'),
  unix(sprintf('mv %s %s; gzip %s', cpf_fname, cpf_fname(1:end-3), cpf_fname(1:end-3))) ;
end ;

fprintf('Done.\n\n') ;
unix(sprintf('rmdir %s', out_dir)) ;
