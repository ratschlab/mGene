function jobinfo =  concat_genomewide_predictons(fn_genome_config,fn_pos,fn_pred,num_splits,avg_all_svms,subtract_diff,run_locally, RPROC)

% concat_genomewide_predictons(fn_genome_config,fn_pos,fn_pred,num_splits,run_locally, RPROC);

jobinfo = rproc_empty(0) ;
q=0;
% options.ncpus = 1;
% options.identifier = sprintf('%sconcat%s_%s',upper(PAR.organism.name(1:4)),PAR.Signal_name,PAR.method.name) ;


RPROC.options.addpaths = {fileparts(which('concat_genomewide_predictons_helper'))};
genome_info = init_genome(fn_genome_config);
P.avg_all_svms = avg_all_svms;
P.subtract_diff = subtract_diff;
P.num_splits = num_splits;
for c = 1:length(genome_info.contig_names)
  for s = '+-'
    P.fn_pred_all = sprintf('%scontig_%i%s_all.mat',fn_pred,c,s);
    P.fn_pred_split = sprintf('%scontig_%i%s_split=',fn_pred,c,s);
    if fexist(P.fn_pred_all)
      continue
    else
      q=q+1;
      if ~run_locally
        [mem_req, time_req, opts] = rproc_memtime_policy('concat_genomewide_predictons_helper',  0, RPROC.options) ;

        if isfield(RPROC, 'collect_jobs') && RPROC.collect_jobs 
          jobinfo(q) = rproc_create('concat_genomewide_predictons_helper', P, mem_req, opts, time_req) ;
        else
          jobinfo(q) = rproc('concat_genomewide_predictons_helper', P, mem_req, opts, time_req);
        end
      else
        concat_genomewide_predictons_helper(P)
      end
    end
  end
end

for c = 1:length(genome_info.contig_names)
  for s = '+-'
    P.fn_pred_all = sprintf('%scontig_%i%s_all.mat',fn_pred,c,s);
    if fexist(P.fn_pred_all)
      vars=who_file(P.fn_pred_all);
      if any(strmatch('Conf',vars)) && any(strmatch('Conf_cum',vars)) && any(strmatch('inventory',vars)) ...
      && any(strmatch('output',vars)) && any(strmatch('pos',vars))
        unix(sprintf('rm -rf %scontig_%i%s_split=*', fn_pred,c,s));
      end
    end
  end
end

if ~run_locally
   rproc_wait(jobinfo, 20, 1, -1);
end
