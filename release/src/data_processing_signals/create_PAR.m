function PAR=create_PAR(Exp, PAR_config, verb)
%  PAR=create_PAR(Exp,PAR_config,verb)
%
% creates PAR with all parameters for one genefinding
% experiment. For a new Organism X a new template file must be added to /template_files/Organisms/X.m. 
% (e.g. C_elegans.m) 
% For each signal a standard method is applied (usually SVM with one or
% more spectrum kernels)
%  
% INPUT  - Exp: filename of current experiment  (e.g. exp_C_elegans_nGASP)
%        - verb: verbose
%
% OUTPUT - PAR. struct with all information needed for genefinding experiment  
  
Signals = {'acceptor','donor','tis','cdsStop','tss','cleave','polya','transacceptor'};
sig_names = {'acc','don','tis','cdsStop','tss','cleave','polya','transacc'};
Methods = {'SVM_WD','SVM_WD','SVM_COMBINED','SVM_COMBINED','SVM_COMBINED','SVM_COMBINED','SVM_COMBINED','SVM_WD'};

content_names = {'intergenic','utr5exon','cds_exon','utr3exon','polya_tail','intron','frame0'};
Methodc = 'SVM_Cont';
changeable_fields = {'organism','label_source','Signals','Contents','LSL','RPROC','SETs','regions','tasks','Source_organisms','FN'} ;

if nargin<3,
  verb=0;
end


PAR = init_PAR;

if nargin>0
  fprintf('parameters from exp file: %s \n',Exp)
  if ~exist('PAR_config', 'var') || isempty(PAR_config),
    PAR_exp = feval(Exp) ;
  else
    try,
      PAR_exp = feval(Exp, PAR_config) ;
    catch
      PAR_exp = feval(Exp) ;
    end 
  end ;
  assert(all(ismember(fieldnames(PAR_exp),changeable_fields)))
  if ~isfield(PAR_exp,'Contents')
    PAR_exp.Contents = struct;
  end
  if ~isfield(PAR_exp,'Signals')
    PAR_exp.Signals = struct;
  end
  PAR.organism = set_default(PAR_exp,'organism',PAR.organism)  ;  
  PAR.label_source = set_default(PAR_exp,'label_source',PAR.label_source) ;
  PAR.RPROC = set_default(PAR_exp,'RPROC',PAR.RPROC)  ; 
  PAR.SETs = set_default(PAR_exp,'SETs',PAR.SETs)  ; 
  PAR.SETs = make_partitions(PAR.SETs);
  PAR.regions = set_default(PAR_exp,'regions',PAR.regions)  ; 
  PAR.tasks = set_default(PAR_exp,'tasks',PAR.tasks)  ; 
  PAR.Source_organisms = set_default(PAR_exp,'Source_organisms',PAR.Source_organisms)  ; 
end

if isempty(PAR.organism.name),
  [organism, label_source] = dummy_organism(PAR.label_source, sig_names);
else
  [organism, label_source] = feval(PAR.organism.name,PAR.label_source,sig_names);
end
PAR.organism = set_default(PAR,'organism',organism)  ;
PAR.label_source = set_default(PAR,'label_source',label_source)  ;

%-----------------------------
%%% prepare all signals
%-----------------------------

for i=1:length(Signals)
  Signal = feval(Signals{i},PAR.label_source);
  if exist('PAR_exp', 'var')
    if isfield(PAR_exp.Signals, sig_names{i}) && isfield(PAR_exp.Signals.(sig_names{i}),'method')
      me = PAR_exp.Signals.(sig_names{i}).method;
      PAR_exp.Signals.(sig_names{i}) = rmfield(PAR_exp.Signals.(sig_names{i}),'method');
    end
    Signal = set_default(PAR_exp.Signals,sig_names{i},Signal)  ;
  end
  Signal.method = feval(Methods{i},Signal);
  if exist('me', 'var')
    PAR_exp.Signals.(sig_names{i}).method = me;
    Signal = set_default(PAR_exp.Signals,sig_names{i},Signal)  ;
    clear me
  end
  % Signal = Signal.domain_adapt.transfer_fct(Signal);
  PAR.Signals.(sig_names{i}) = Signal;
end

%-----------------------------
%%% prepare all contents
%-----------------------------
for i=1:length(content_names)
  Content = feval(content_names{i});
  if exist('PAR_exp', 'var')
    if isfield(PAR_exp.Contents, content_names{i}) && isfield(PAR_exp.Contents.(content_names{i}),'method')
      me = PAR_exp.Contents.(content_names{i}).method;
      PAR_exp.Contents.(content_names{i}) = rmfield(PAR_exp.Contents.(content_names{i}),'method');
    end
    Content = set_default(PAR_exp.Contents,content_names{i},Content)  ;
  end
  Content.method = feval(Methodc,Content);
  if exist('me', 'var')
    PAR_exp.Contents.(content_names{i}).method = me;
    Content = set_default(PAR_exp.Contents,content_names{i},Content)  ;
    clear me
  end
  % Content = Content.domain_adapt.transfer_fct(Content);
  PAR.Contents.(content_names{i}) = Content;
end

%-----------------------------
%%% prepare LSL
%-----------------------------

LSL = lsl_template; 
% LSL = LSL.domain_adapt.transfer_fct(LSL);
if exist('PAR_exp', 'var')
  LSL = set_default(PAR_exp,'LSL',LSL)  ; 
end
PAR.LSL=LSL;

if isfield(PAR_config, 'track_names')&&isfield(PAR_config, 'segment_feature_names')
  PAR.LSL.method.track_names 				= PAR_config.track_names;
  PAR.LSL.method.track_functions 			= PAR_config.track_functions;
  PAR.LSL.method.track_files 				= PAR_config.track_files;
  PAR.LSL.method.track_params 				= PAR_config.track_params;
  PAR.LSL.method.track_monoton_functions 		= PAR_config.track_monoton_functions;

  PAR.LSL.method.segment_feature_names 			= PAR_config.segment_feature_names;
  PAR.LSL.method.segment_feature_functions 		= PAR_config.segment_feature_functions;
  PAR.LSL.method.segment_feature_files 			= PAR_config.segment_feature_files;
  PAR.LSL.method.segment_feature_params 		= PAR_config.segment_feature_params;
  PAR.LSL.method.segment_feature_monoton_functions 	= PAR_config.segment_feature_monoton_functions;
else
  PAR.LSL.method.track_names 				= {}; 
  PAR.LSL.method.track_functions 			= {}; 
  PAR.LSL.method.track_files 				= {}; 
  PAR.LSL.method.track_params 				= {}; 
  PAR.LSL.method.track_monoton_functions 		= {}; 

  PAR.LSL.method.segment_feature_names 			= {}; 
  PAR.LSL.method.segment_feature_functions 		= {}; 
  PAR.LSL.method.segment_feature_files 			= {}; 
  PAR.LSL.method.segment_feature_params 		= {}; 
  PAR.LSL.method.segment_feature_monoton_functions 	= {}; 

end

%-----------------------------
%%% prepare filenames
%-----------------------------
FN = init_FN;
if  exist('PAR_exp', 'var')
  FN = set_default(PAR_exp,'FN',FN)  ;
end

FN = create_filenames(FN,PAR.organism,PAR.Signals,PAR.Contents,verb);
if  exist('PAR_exp', 'var')
  FN = set_default(PAR_exp,'FN',FN)  ;
end
PAR.FN=FN;

%-----------------------------
%%% Add polya consensus if file available
%-----------------------------

if fexist(PAR.FN.input_sig.polya.fn_polya)
  load(PAR.FN.input_sig.polya.fn_polya,'consensus')
  PAR.Signals.polya.consensus = consensus;
end

%PAR.FN.source_directory = '~/svn/projects/genefinding/';
%[tmp,PAR.FN.svn_info]=unix('svn info');

PAR.FN.svn_info = '' ;

fprintf('loading model from %s\n', PAR.LSL.method.model_fct)
PAR.model = feval(PAR.LSL.method.model_fct,PAR.organism,PAR.LSL, ...
                  PAR.Signals);

