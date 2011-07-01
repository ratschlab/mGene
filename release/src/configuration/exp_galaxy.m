function PAR=exp_galaxy(PAR_config)
% PAR=exp_galaxy(PAR_config)

PAR.SETs.num_partitions=5 ;
PAR.SETs.num_splits=5 ;

PAR.FN.input_dir=PAR_config.dir ;
PAR.FN.output_dir=PAR_config.dir;
PAR.FN.genome_dir=PAR_config.dir ;
PAR.FN.input.fn_genome_config = PAR_config.genome_config ;

% preparing RPROC settings
PAR.RPROC.MEMREQ = 6000 ;
PAR.RPROC.collect_jobs = 0 ;
PAR.RPROC.exm_per_batch = 1 ;
PAR.RPROC.options.addpaths = {} ;

signals = {'acc','don','tis','cdsStop','tss','cleave'};
for i=1:length(signals)
  PAR.label_source.(signals{i}).from_fn_candsites = 1 ;
  PAR.Signals.(signals{i}).filter_label.train.USE_ALL = 1;
  PAR.Signals.(signals{i}).filter_label.test.USE_ALL = 1;
  PAR.Signals.(signals{i}).filter_label.eval.USE_ALL = 1;
  PAR.Signals.(signals{i}).Conf_names = {} ;
  PAR.Signals.(signals{i}).method.threads = 1;
end


if 0,
  PAR.Signals.acc.order = 10 ;
  PAR.Signals.acc.shift = 0 ;

  PAR.Signals.don.order = 10 ;
  PAR.Signals.don.shift = 0 ;
  
  PAR.Signals.tis.order =  {[],[5],[]};
  PAR.Signals.tis.shift =  {[],[0],[]};
  
  PAR.Signals.cdsStop.order = {[],[5],[]};
  PAR.Signals.cdsStop.shift = {[],[0],[]};
  
  PAR.Signals.tss.order = {[],[5],[]};
  PAR.Signals.tss.shift = {[],[0],[]};
  
  PAR.Signals.cleave.order = {[],[5],[]};
  PAR.Signals.cleave.shift = {[],[0],[]};

end


contents = {'intergenic','utr5exon','cds_exon','utr3exon','polya_tail','intron','frame0'};

for i=1:length(contents)
  PAR.label_source.(contents{i}).from_fn_candsites = 1 ;
  PAR.Contents.(contents{i}).filter_label.train.USE_ALL = 1;
  PAR.Contents.(contents{i}).filter_label.test.USE_ALL = 1;
  PAR.Contents.(contents{i}).filter_label.eval.USE_ALL = 1;
  PAR.Contents.(contents{i}).Conf_names = {} ;
  PAR.Contents.(contents{i}).method.threads = 1;
end

PAR.LSL.method.use.signals.polya = 0 ;
PAR.LSL.method.exm_per_solve = 100 ;

PAR.regions.max_length = 30000 ;


if 1,

%----- REGULARIZATION
method=PAR.LSL.method ;

%C = 1 ;
C = 1e-3;
%%% linear regularization
method.par_ms.C_regul.lengths = C ;        % linear smoothness for len_hist
method.par_ms.C_regul.signals = C ;        % for sig_hist
method.par_ms.C_regul.contents = C ;       % for content_hist
method.par_ms.C_regul.transitions = C ;    % for

method.par_ms.C_regul.tiling = C ;  
method.par_ms.C_regul.rna_seq = C ;  


%%% quadratic regularization
C_sq = 1e-3;

method.par_ms.C_regul.transitions_sq = C_sq;
method.par_ms.C_regul.plif_ys_sq = C_sq*0.1 ;  % small_plif_ys_penalty
method.par_ms.C_regul.smoothness_sq = C_sq ;

PAR.LSL.method=method ;

end
