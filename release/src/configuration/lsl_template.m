function LSL = lsl_template
  

method.name = 'lsl' ;
method.model_fct = 'make_model';
method.train_fct = 'lsl_train';
method.predict_fct = 'lsl_predict';

method.dt = datestr(now,'dd-mm-yyyy');
method.user = 'GS';
method.plif_bins   = 20 ;   % number of plif bins

%%%

%%Signals
use.signals.trans = 0 ;
use.signals.transacc = 0 ;
use.signals.tss = 1;
use.signals.polya = 1;
use.signals.cleave =  1;

%%States
use.states.SL1 = 0 ;
use.states.SL2 = 0;

if 0
	warning('use non-coding model')
	use.non_coding = 1;
end

%%transitions
use.transitions.allow_missing_tss = 0;
use.transitions.one_gene = 0 ; % force exactly one gene
use.transitions.use_multi = 1 ; % allow zero or more than one gene 

%%loss
use.loss.orf_only = 0;
use.model.one_gene = 0;
use.model.allow_missing_tss = 0;
use.model.use_multi = 1;
%use.nucleotide_loss_factor = 1e-2;

%%Segments & lengths
use.segments.transexon = 0 ;
use.segments.intercistronic = 0 ;
use.segments.intergenictrans = 0;
use.segments.polya_tail = 1;

%%Content 

use.contents.frame = 1;
use.contents.pre_comp = 1;


%%rna_seq_polya
if 0
	use.segments.rna_seq_polya = 1;
	use.contents.rna_seq_polya = 1;
end


%%additional features
use.addfeatures.conservation = 0 ;
use.addfeatures.ests = 0 ;
use.addfeatures.proteins = 0 ;
%use.addfeatures.tiling = 0;
%use.addfeatures.rna_seq = 0;
%use.addfeatures.rna_seq_intron_track = 0;
%use.addfeatures.rna_seq_intron_list = 0;
%use.addfeatures.rna_seq_intron_quality = 0;

% keyboard
method.use = use;



%----- REGULARIZATION
%C = 1 ;
%C = 1e-1;
C = 0;
%%% linear regularization
method.par_ms.C_regul.lengths = C ;        % linear smoothness for len_hist
method.par_ms.C_regul.signals = C ;        % for sig_hist
method.par_ms.C_regul.contents = C ;       % for content_hist
method.par_ms.C_regul.transitions = C ;    % for

method.par_ms.C_regul.tiling = C ;  
method.par_ms.C_regul.rna_seq = C ;  
method.par_ms.C_regul.rna_seq_intron_track = C ;  
method.par_ms.C_regul.rna_seq_intron_list = C ;  
method.par_ms.C_regul.rna_seq_intron_quality = C ;  


%%% quadratic regularization
C_sq = 1e-1;

method.par_ms.C_regul.transitions_sq = C_sq;
method.par_ms.C_regul.plif_ys_sq = C_sq*0.1 ;  % small_plif_ys_penalty
method.par_ms.C_regul.smoothness_sq = C_sq ;


%----------- for predictions only
method.par_ms.cleave_offset = 6;
%-----------

method.solver = 'dual';
method.ignore_init_file = 1;

method.opt_step = 1000 ;
method.max_neg = 2 ;
method.margin = 1 ;

method.INF = 1e20;
method.max_num_iterations = 100;

method.exm_per_solve = 400 ;
% method.exm_per_batch = 5 ;

%method.cleave_offset = 6 ;
% -------
% -------
LSL.method = method;

