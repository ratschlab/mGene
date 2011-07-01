function jobinfo = convert_predictions2confs(fn_genome_config,fn_pred,fn_svms,num_splits,run_locally, RPROC)
% jobinfo = convert_predictions2confs(fn_genome_config,fn_pred,fn_svms,num_splits,run_locally, RPROC)

jobinfo = rproc_empty(0) ;
% options.ncpus = 1;
% options.identifier = sprintf('%sconf%s_%s',upper(PAR.organism.name(1:4)),PAR.Signal_name,PAR.method.name) ;

fid=1;q=0;
for p = 1:num_splits
  fprintf(fid,'partition %i\n',p);
  
  fn_SVM = sprintf('%s%i.mat',fn_svms,p) ;
 
  %--------------------------------------------------------------------------------%
  % load parameters from SVM file
  %--------------------------------------------------------------------------------%
  
  genome_info = init_genome(fn_genome_config);
  clear  prob prob_cum limits Conf Conf_cum

  if ~(all(ismember({'prob', 'prob_cum','limits'}, who_file(fn_SVM))))
    warning('Variables prob, prob_cum and or limits not found in %s',fn_SVM)
    error('calculate PLIFS first: execute this function with PAR.tasks.pred_from_seq=0; to calculate PLIFS')
  else
    load(fn_SVM, 'Train');
    PAR_new = load(fn_SVM, 'prob', 'prob_cum', 'limits');
    if ~isfield(Train.PAR.SETs, 'which_pred')
      Train.PAR.SETs.which_pred = Train.PAR.SETs.which_pred_svm;
    end
  end
  assert(Train.PAR.SETs.which_pred == p)
  %------------------------------------------------------------------------------------%
  % perform transformations
  %------------------------------------------------------------------------------------%
  for c = 1:length(genome_info.contig_names)
    for s = '+-'
      PAR_new.fn_pred = sprintf('%scontig_%i%s_split=%i.mat',fn_pred,c,s,p);
      for j=1:20
        if ~fexist(PAR_new.fn_pred)
	  fprintf(fid,'\n waiting for file: %s\n',PAR_new.fn_pred);
          pause(30) 
	end
      end
      if ~fexist(PAR_new.fn_pred)
       fprintf(fid,'\n file not found: %s', PAR_new.fn_pred);
       fprintf(fid,'\n run genomewide_prediction first\n');
      %elseif all(ismember({'Conf', 'Conf_cum'}, who_file(PAR_new.fn_pred))),
      %  continue
      else
        fprintf(fid,'\n start transforming on contig %i%s\n',c,s);
        %%% Convert SVM Outputs on test set `
        q=q+1;
        if ~run_locally

          [mem_req, time_req, opts] = rproc_memtime_policy('Out2Confhelper',  0, RPROC.options) ;

          if isfield(RPROC, 'collect_jobs') && RPROC.collect_jobs 
            jobinfo(q) = rproc_create('Out2Confhelper', PAR_new, mem_req, opts, time_req) ;
          else
            jobinfo(q) = rproc('Out2Confhelper', PAR_new, mem_req, opts, time_req);
          end
        else
          Out2Confhelper(PAR_new);
        end
      end
    end%loop over strands
  end%loop over contigs
end %% loop over partitions

rproc_wait(jobinfo, 20, 1, -1);
