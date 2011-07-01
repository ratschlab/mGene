function PAR=init_PAR(signal_names,content_names);
 
  
% -------------------------------------------------------------------------------------
% Experiments
% -------------------------------------------------------------------------------------
if nargin<1
  signal_names = {'acc','don','tis','cdsStop','tss','cleave','polya','transacc'};
end
if nargin<2
  content_names = {'intergenic','utr5exon','cds_exon','utr3exon','polya_tail','intron','frame0'};
end




PAR.date = date;
PAR.user = whoami;
  
PAR.organism.clade = [];
PAR.organism.name = '';
PAR.organism.full_name = [];
PAR.organism.release = [];
PAR.organism.version = [];
PAR.organism.genebuilt = [];

% -------------------------------------------------------------------------------------
% info regarding original gene structure
% -------------------------------------------------------------------------------------


info_genes.transcript.anno = -1;
info_genes.transcript.est = 0;
info_genes.transcript.cDNA = 1;
info_genes.transcript.fulllength = 2;

info_genes.transcript_status.confirmed = 1;
info_genes.transcript_status.partially_confirmed = 0;
info_genes.transcript_status.predicted = -1;


info_genes.pseudo.dummy = [];

% -------------------------------------------------------------------------------------
% info regarding label source
% -------------------------------------------------------------------------------------

label_source.use_anno = [];
label_source.use_confirmed = [];
for j=1:length(signal_names)
  signal_name = signal_names{j};
  label_source.(signal_name).anno = [];
  label_source.(signal_name).est = [];
  label_source.(signal_name).cDNA = [];
  label_source.(signal_name).fulllength = [];
  label_source.(signal_name).db = [];
  label_source.(signal_name).from_fn_candsites = [];
  
  info_genes.(signal_name).anno = -1;
  info_genes.(signal_name).est = 0;
  info_genes.(signal_name).cDNA = 1;
  info_genes.(signal_name).fulllength = 2;
  info_genes.(signal_name).db = 3;
  info_genes.(signal_name).from_fn_candsites = 4;
end
for j=1:length(content_names)
  content_name = content_names{j};
  label_source.(content_name).from_fn_candsites = [];
end

label_source.tis.maxorf = [];
info_genes.tis.maxorf = 5;
label_source.cdsStop.maxorf = [];
info_genes.cdsStop.maxorf = 5;

label_source.transacc.anno_SL1 = [];
label_source.transacc.anno_SL2 = [];
label_source.transacc.pred = [];
info_genes.transacc.anno_SL1 = 5;
info_genes.transacc.anno_SL2 = 6;
info_genes.transacc.pred = 7;

label_source.polya.from_cleave = [];
info_genes.polya.from_cleave = 5;


PAR.label_source = label_source;
PAR.info_genes = info_genes;


PAR.Signal_name = [];
for j=1:length(signal_names)
  signal_name = signal_names{j};
  PAR.Signals.(signal_name) = struct;
end

cnt=0;
PAR.Content_name = [];
for j=1:length(content_names)
  PAR.Contents.(content_names{j}) = struct;
  PAR.Contents.(content_names{j}).name = content_names{j};
  PAR.Contents.(content_names{j}).type = cnt;
  cnt=cnt+1;
end

PAR.LSL.method = struct;

% -------------------------------------------------------------------------------------
% memory requirements
% -------------------------------------------------------------------------------------

PAR.RPROC.collect_jobs = 0 ; 
PAR.RPROC.options.ncpus = 1 ; 
PAR.RPROC.options.waitonfull = [];
PAR.RPROC.options.maxjobs = [] ;
PAR.RPROC.options.immediately = 0;
PAR.RPROC.options.immediately_bg = 0;
PAR.RPROC.options.addpaths = [];
PAR.RPROC.options.rmpaths = [];

PAR.RPROC.MEMREQ = [];
PAR.RPROC.time_req = [];
PAR.RPROC.exm_per_batch = 1 ;


% -------------------------------------------------------------------------------------
% SPLITS AND PARTITIONS
% -------------------------------------------------------------------------------------

PAR.SETs.num_splits = 5;
PAR.SETs.num_partitions = 5;
PAR.SETs.train_subsample_neg = 1;
PAR.SETs.test_subsample_neg = 1;
PAR.SETs.train_subsample_pos = 1;
PAR.SETs.test_subsample_pos = 1;
%for lsl
PAR.SETs.num_train = [] ; 
PAR.SETs.ignore_incomplete = 1 ; 


% -------------------------------------------------------------------------------------
% original data processing stuff
% -------------------------------------------------------------------------------------

PAR.regions.max_length = 1000000;
PAR.regions.min_genes_per_region = 1;
PAR.regions.offset_test = 20000;
PAR.regions.offset_train = 1000;
PAR.regions.minCutDist = 250;
PAR.regions.cluster_paralogs = 0;
PAR.regions.merge_confirmed = 0;

% -------------------------------------------------------------------------------------
% TASKS
% -------------------------------------------------------------------------------------

PAR.tasks.signals.prepare_labels = 0;
PAR.tasks.signals.train = 0;
PAR.tasks.signals.modelsel = 0;
PAR.tasks.signals.learn_plifs = 0;
PAR.tasks.signals.get_train_error = 0;

PAR.tasks.signals.pred_on_genome = 0; 
PAR.tasks.signals.verbose = 0; 

PAR.tasks.signals.run_locally = 0;

% -------------------------------------------------------------------------------------
% transfer from other organisms
% -------------------------------------------------------------------------------------

PAR.Source_organisms.names = {};
PAR.Source_organisms.dir_names = {};
PAR.Source_organisms.exp_names = {};


% -------------------------------------------------------------------------------------
% FN
% -------------------------------------------------------------------------------------
PAR.FN = struct;

