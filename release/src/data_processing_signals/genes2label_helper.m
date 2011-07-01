function genes2label_helper(genes, genome_info, genome_config, spf_fname, spf_dir, signal_name, covered_genic_positions_only) ;
% genes2label_helper(genes, genome_info, genome_config, spf_fname, spf_dir, signal_name, covered_genic_positions_only) ;

disp('-------------------------------------------------------') ;
disp('Step 5: Extracting labeled positions from the genome...') ;
disp('-------------------------------------------------------') ;

if ~isempty(spf_dir),
  galaxy_PAR.dir = spf_dir ;
else
  galaxy_PAR.dir = tempname ;
end ;
galaxy_PAR.genome_config = genome_config ;
PAR = create_PAR('exp_galaxy', galaxy_PAR);
PAR.Signal_name = signal_name ;
PAR.Signals.(signal_name).Conf_names = {} ;


if ~isempty(spf_dir),
  out_dir = [spf_dir '/tmp'] ;
else
  out_dir = tempname ;
end ;
unix(sprintf('mkdir -p %s', out_dir)) ;

contig_names = genome_info.contig_names ;
for c=1:length(contig_names)
  for s='+-',
    cs_idx = find(ismember(upper({genes.chr}), upper(contig_names{c})) & [genes.strand]==s) ;
    cs_genes = genes(cs_idx) ;

    if covered_genic_positions_only,
      covered_genic_positions = {} ;
      for i=1:length(cs_genes),
        covered_genic_positions{end+1} = cs_genes(i).start-1:cs_genes(i).stop+1 ;
      end ;
      covered_genic_positions=[covered_genic_positions{:}] ;
    end ;

    % get the list of positively labeled positions (those in the annotation)
    seq = load_genomic(contig_names{c}, s, 1, inf, genome_info) ;

    P_=struct ;
    P_.genes = cs_genes;
    P_.seq = seq;
    P_.Signal = PAR.Signals.(PAR.Signal_name);
    P_.info = PAR.info_genes.(PAR.Signal_name);
    P_.signal_name = PAR.Signal_name ;
    P_.strand = s; 
    [pos, label] = feval(PAR.Signals.(signal_name).label_fct, P_);
    clear P_ ;

    if covered_genic_positions_only,
      [tmp,idx1,idx2]=intersect(covered_genic_positions, pos) ;
      pos   = pos(idx2) ;
      label = label(idx2) ;
    end ;

    
    %subsample negative examples
    pos_idx = find(label==1);
    neg_idx = find(label==-1);
    if length(pos_idx)>0 && length(neg_idx)>0
      real_ratio = length(pos_idx)/length(neg_idx) ;
      fprintf('current class ratio (pos/neg): %2.1f%%\n', real_ratio) ;
      if (PAR.Signals.(signal_name).sampling.pos_neg_ratio - real_ratio) > 0
        % if real_ratio is smaller than pos_neg_ratio, there are to many negative examples
        rand('state', 24593854);
        keep = randperm(length(neg_idx));
        neg_idx = neg_idx(keep(1:floor(length(pos_idx) / PAR.Signals.(signal_name).sampling.pos_neg_ratio)));
        pos = pos([pos_idx neg_idx]);
        label = label([pos_idx neg_idx]);
        real_ratio = length(pos_idx)/length(neg_idx);
        fprintf('class ratio (pos/neg) after downsampling: %2.1f%%\n', real_ratio) ;
      end 
    end

    % sort the positions
    [tmp,idx] = sort(pos) ;
    pos = pos(idx) ;
    label = label(idx) ;
    
    % save them to a file
    filename_prefix = sprintf('%s/contig_%i%c', out_dir, c, s) ;
    fn_in=[filename_prefix '_label.spf'];
    unix(sprintf('rm -f %s', fn_in)) ;
    save_predictions_as_sigpred(filename_prefix, pos, label, 'label', signal_name, contig_names{c}, s, 0) ;
    save_predictions_as_sigpred_bin(filename_prefix, pos, label, 'label', signal_name, contig_names{c}, s) ;
  end ;
end ;
fprintf('Done.\n\n') ;

disp('----------------------------------') ;
disp('Step 6: Generating output files...') ;
disp('----------------------------------') ;
unix(sprintf('rm -f %s', spf_fname)) ;
cnt = 0 ;
for c = 1:length(contig_names)
  for s = '+-'
    P.chrom = c;
    P.strand = s;
    filename = sprintf('%s/contig_%i%s', out_dir, c, s);
    fn_in = [filename '_label.spf'];
    %ret = unix(['cat ' fn_in ' >> ' spf_fname '; rm -f ' fn_in]) ;
    %assert(ret==0) ;

    fn_mat = sprintf('%s/contig_%i%s_label_spf.mat', out_dir, c, s);
    fnall_mat = sprintf('%s/label_spf.mat', spf_dir);
    SPF=load(fn_mat);
    %fieldnames(SPF)

    save_append(fnall_mat, 1, sprintf('SPF_%i', cnt), SPF) ;
    cnt = cnt + 1 ;
    ret = unix(['rm -f ' fn_mat]) ;
    assert(ret==0) ;
  end ;
end ;
if length(spf_fname)>=3 && isequal(spf_fname(end-2:end), '.gz'),
  unix(sprintf('mv %s %s; gzip %s', spf_fname, spf_fname(1:end-3), spf_fname(1:end-3))) ;
end ;

fprintf('removing directory %s\n', out_dir);
unix(sprintf('rm -r %s', out_dir)) ;
fprintf('Done.\n\n') ;
