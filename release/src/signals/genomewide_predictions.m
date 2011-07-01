function jobinfo = genomewide_predictions(fn_genome_config,fn_pos,fn_pred,fn_svms,num_splits,avg_all_svms,run_locally, RPROC, verb)
% jobinfo = genomewide_predictions(fn_pos,fn_pred,fn_genome_config,verb)


if nargin<9
  verb=0;
end
if nargin<7
  run_locally=1;
end
if nargin<6
  avg_all_svms = 0;
end

jobinfo = rproc_empty(0) ;
if ~run_locally
  jobinfo = rproc_empty(0) ;
  options = RPROC.options;
  options.ncpus = 1;
  [sg_path, envstr] = shogun_settings();
  options.start_dir = get_base_dir();
  options.addpaths = {sg_path, fileparts(which('gen_all_positions_rproc'))};
  options.envstr = envstr;

  % options.identifier = sprintf('%spred%s_%s',upper(PAR.organism.name(1:4)),PAR.Signal_name,PAR.method.name) ;
  MEMREQ = RPROC.MEMREQ;
end

    
fid=1;
q=0;

genome_info = init_genome(fn_genome_config);
num_pred = num_splits * length(genome_info.contig_names)*2;
fprintf(fid,'total number of jobs for prediction (num_svms*num_contigs*num_strands): %i\n',num_pred);

trivial_regions = init_regions(fn_genome_config);

for p = 1:num_splits
  fprintf(fid,'partition %i\n',p);
  
  fn_SVM = sprintf('%s%i.mat',fn_svms,p) ;
  %--------------------------------------------------------------------------------%
  % load parameters from SVM file
  %--------------------------------------------------------------------------------%
  
  load(fn_SVM,'Train');
  assert(Train.PAR.SETs.which_pred == p)
  
  PAR_new = Train.PAR;
  
  if isfield(PAR_new.Signal.method,'genomewide_predict_fct') && ~isempty(PAR_new.Signal.method.genomewide_predict_fct)
    PAR_new.Signal.method.predict_fct =  PAR_new.Signal.method.genomewide_predict_fct;
  end
  PAR_new.fn_SVMs = fn_SVM;
  PAR_new.pred_from_seq=1;
  PAR_new.avg_all_svms = avg_all_svms ; 
  if ~isfield(PAR_new.SETs,'which_pred')
    PAR_new.SETs.which_pred=PAR_new.SETs.which_pred_svm;
  end  
  assert(Train.PAR.SETs.which_pred == p)
  clear Train
  
  %------------------------------------------------------------------------------------%
  % perform genome-wide predictions
  %------------------------------------------------------------------------------------%
  for r_idx = 1:length(trivial_regions) 
    c = trivial_regions(r_idx).chr_num;
    s = trivial_regions(r_idx).strand;
    PAR_new.fn_pos = fn_pos;
    PAR_new.fn_pred = sprintf('%scontig_%i%s_split=%i.mat',fn_pred,c,s,p);
    if verb
      fprintf(fid,'\n start predicting on contig %i%s\n',c,s);
    end
    if ~fexist(PAR_new.fn_pred) 
      PAR_new.region = trivial_regions(r_idx);
      PAR_new.plif_on = 0;
      PAR_new.eval_on = 0;

      q=q+1; 
      if ~run_locally


        if ~isfield(RPROC,'divide4pred') || ~RPROC.divide4pred

          %[mem_req, time_req, opts] = rproc_memtime_policy('load_and_predict',  0, options) ;
          [mem_req, time_req, opts] = rproc_memtime_policy(sprintf('load_and_predict:%s', PAR_new.Signal.name),  0, options) ;

          if isfield(RPROC, 'collect_jobs') && RPROC.collect_jobs 
            jobinfo(q) = rproc_create('load_and_predict', PAR_new, mem_req, opts, time_req) ;
          else
            jobinfo(q) = rproc('load_and_predict', PAR_new, mem_req, opts, time_req) ;
          end
        else
          q=q-1;
          
          output = [];
          pos = [];
          all_done=1;Temp_old =[];
          for chunk_idx=1:RPROC.numchunks
            PAR_new.chunk_idx = chunk_idx;PAR_new.RPROC.numchunks = RPROC.numchunks ;
            PAR_new.RPROC.divide4pred = RPROC.divide4pred;

            [mem_req, time_req, opts] = rproc_memtime_policy('load_and_predict',  0, options) ;

            if  ~fexist(sprintf('%s%i.mat',PAR_new.fn_pred, chunk_idx)) && isfield(RPROC, 'collect_jobs') && RPROC.collect_jobs 
              q=q+1;
              jobinfo(q) = rproc_create('load_and_predict', PAR_new, mem_req, opts, time_req);
              all_done=0;
            elseif  ~fexist(sprintf('%s%i.mat',PAR_new.fn_pred, chunk_idx))
              q=q+1;
              jobinfo(q) = rproc('load_and_predict', PAR_new, mem_req, opts, time_req);
              all_done=0;
            else 
              Temp = load(sprintf('%s%i.mat',PAR_new.fn_pred, chunk_idx) )
              Temp.PAR = rmfield(Temp.PAR,'RPROC');
              Temp.PAR = rmfield(Temp.PAR,'chunk_idx');
              if ~isempty(Temp_old)
                assert(isequal(Temp_old.svmFileName,Temp.svmFileName));
                assert(isequal(Temp_old.PAR,Temp.PAR));
                assert(isequal(Temp_old.sg_ver,Temp.sg_ver));
              end
              Temp_old = Temp;
              output = [output Temp.output];
              pos = [pos Temp.pos];
            end
          end
          if all_done
            Temp = rmfield(Temp,'output');Temp = rmfield(Temp,'pos');
            %save(PAR_new.fn_pred,'struct','Temp');
            save_append(PAR_new.fn_pred, 0, Temp);
            %save(PAR_new.fn_pred,'-append','output','pos');
            inventory = {'output', 'pos'};
	    assert(isequal(size(pos),size(output)))
            save_append(PAR_new.fn_pred, 1, 'inventory', inventory,  'output', output, 'pos', pos);
          end
        end
      else
        load_and_predict(PAR_new)
      end
    else
      if verb
        fprintf(fid,'\n file already exists ');
      end
    end
  end%loop over contigs
end %% loop over partitions

if ~exist('jobinfo', 'var')| isempty(jobinfo),
  jobinfo = rproc_empty(0) ;
  return 
end

% [jobinfo, meta_jobinfo] = rproc_submit_batch(jobinfo, PAR.RPROC.exm_per_batch) ;
if ~(isfield(RPROC, 'collect_jobs') && RPROC.collect_jobs),
%   [jobinfo,num_crashed] = rproc_wait(jobinfo,20,1,0);
end ;

[jobinfo,num_crashed] = rproc_wait(jobinfo, 30, 1, -1);

