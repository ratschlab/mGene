function [jobinfo] = save_predictions(fn_genome_config, fn_pred, ...
                                      Signal_name, resolution, conf_cum_thresh, strands, run_locally, RPROC, do_gzip, save_as, inventory)
%  [jobinfo] =
%  save_predictions(fn_genome_config,fn_pred,strands,run_locally,RPROC,do_gzip,save_as,inventory)
%
%  Loads predictions corresponding to each contig for a signal
%  defined in PAR.Signal_name.
%  Saves them to binary files which can be read by 
%  interval_query.
%

%
%  see retrieve_signals, interval_query 
if nargin<2
  strands='+-';
end

if ~exist('do_gzip', 'var'),
  do_gzip = 1 ;
end ;
if ~exist('save_as', 'var'),
  save_as.spf_ascii = 1 ;
  save_as.spf_binary = 1 ;
  save_as.wiggle = 0 ;
end ;

jobinfo = rproc_empty(0) ;
% options = RPROC.options;
% options.ncpus = 1;
% options.identifier = sprintf('%scheck%s_%s',upper(PAR.organism.name(1:4)),PAR.Signal_name,PAR.method.name) ;
% MEMREQ = RPROC.MEMREQ;

fid=1;q=0;

genome_info = init_genome(fn_genome_config);

P.fn_genome_config = fn_genome_config;
P.Signal_name = Signal_name;
P.fn_pred = fn_pred;
P.conf_cum_thresh = conf_cum_thresh;
P.resolution = resolution ;
P.save_as = save_as ;


% options.identifier = sprintf('%ssave%s_%s',upper(PAR.organism.name(1:4)),PAR.Signal_name,PAR.method.name) ;
for c = 1:length(genome_info.contig_names)
  for s = strands
    P.chrom = c;
    P.strand = s;
    P.do_gzip = do_gzip ;
    filename = sprintf('%scontig_%i%s', fn_pred, c, s);

    one_missing=0 ;
    for k=1:length(inventory)
      if ~fexist([filename '.' inventory{k}]) || ~fexist([filename '_' inventory{k} '.spf']) || ~fexist([filename '_' inventory{k} '_spf.mat'])
        one_missing = 1 ;
      end ;
    end ;
    if one_missing,
      if ~run_locally
        q=q+1;

        [mem_req, time_req, opts] = rproc_memtime_policy('save_predictions_helper',  0, RPROC.options) ;%RPROC.MEMREQ ;

        if isfield(RPROC, 'collect_jobs') && RPROC.collect_jobs 
          jobinfo(q) = rproc_create('save_predictions_helper', P, mem_req, opts, time_req) ;
        else
          jobinfo(q) = rproc('save_predictions_helper', P, mem_req, opts, time_req);
        end
      else
        save_predictions_helper(P)
      end
    end
  end %%% loop over strand
end %%% loop over contigs 


rproc_wait(jobinfo, 20, 1, -1);


