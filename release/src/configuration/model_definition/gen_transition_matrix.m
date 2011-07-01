function model = gen_transition_matrix(model,PAR)

use=model.use;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define transitions

A = -inf(model.cnt_states, model.cnt_states);
p = -inf(model.cnt_states,1);
q = -inf(model.cnt_states,1) ;
p(model.state_ids.seq_start) = 0 ;
q(model.state_ids.seq_end) = 0 ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start constructing transition matrix A

orf_from = -1*ones(model.cnt_states,1) ;
orf_to = -1*ones(model.cnt_states,1) ;

orf_from(model.state_ids.acc(2:7))= [ 0 1 1 2 2 2] ;
orf_from(model.state_ids.tis)     = 0 ; % AGAGT case
if model.use.states.SL1 & model.use.states.SL2
  orf_from(model.state_ids.SL1(2)) = 0 ; % AGAGT case
  orf_from(model.state_ids.SL2(2)) = 0 ; % AGAGT case
elseif model.use.signals.trans
  orf_from(model.state_ids.trans(2)) = 0 ; % AGAGT case
end
orf_to(model.state_ids.don(2:7))  = [0 1 1 2 2 2] ;
orf_to(model.state_ids.cdsStop(1))   = 0 ;
orf_info=[orf_from, orf_to] ;

%%%%%%%%%%%
%% no gene / intergenic stuff
%%%%%%%%%%%
if ~use.model.one_gene,
  A(model.state_ids.seq_start,  model.state_ids.seq_end) = model.lengths.intergenic ;
end
%A(model.state_ids.intergenic_long, model.state_ids.intergenic_long) = 0; % no penalty et al
%A(model.state_ids.seq_start, model.state_ids.intergenic_long) = model.lengths.intergenic ;
%A(model.state_ids.intergenic_long, model.state_ids.seq_end) = model.lengths.intergenic ;

%%%%%%%%%%%
% 5' UTR 
%%%%%%%%%%%
gene_start_states_in = [];
gene_start_states_out = [];
mRNA_start_states_out = [model.state_ids.tis];

if use.signals.tss,
  gene_start_states_in = [gene_start_states_in model.state_ids.tss] ;
  gene_start_states_out = [gene_start_states_out model.state_ids.tss] ;
  if use.model.allow_missing_tss,
    assert(false) ;
    gene_start_states_in = [gene_start_states_in model.state_ids.tis] ;
  end ;
else
  assert(false)
  gene_start_states_in = [gene_start_states_in model.state_ids.tis] ;
end ;

if isfield(model.use, 'non_coding')&&model.use.non_coding
	A(model.state_ids.tss_nc, model.state_ids.don_nc) = model.lengths.first_nc_exon;
	A(model.state_ids.acc_nc, model.state_ids.cleave_nc) = model.lengths.last_nc_exon;
	A(model.state_ids.don_nc, model.state_ids.acc_nc) = model.lengths.intron;
	A(model.state_ids.acc_nc, model.state_ids.don_nc) = model.lengths.middle_nc_exon;
	A(model.state_ids.cleave_nc, model.state_ids.tss) = model.lengths.intergenic;
	A(model.state_ids.cleave_nc, model.state_ids.tss_nc) = model.lengths.intergenic;
	A(model.state_ids.cleave, model.state_ids.tss_nc) = model.lengths.intergenic;
	A(model.state_ids.seq_start, model.state_ids.tss_nc) = model.lengths.intergenic;
	A(model.state_ids.cleave_nc, model.state_ids.seq_end) = model.lengths.intergenic;
	A(model.state_ids.tss_nc, model.state_ids.cleave_nc) = model.lengths.first_nc_exon;
end

if use.signals.trans,
  gene_start_states_in = [gene_start_states_in model.state_ids.SL1] ;
  gene_start_states_out = [gene_start_states_out model.state_ids.SL1(1)] ;
  gene_start_states_out = [gene_start_states_out model.state_ids.SL2(1)] ;
  mRNA_start_states_out = [mRNA_start_states_out model.state_ids.SL1(2)] ;%'AGAGT cases'
  mRNA_start_states_out = [mRNA_start_states_out model.state_ids.SL2(2)] ;%'AGAGT cases'
  if use.signals.tss, 
    A(model.state_ids.tss, model.state_ids.SL1)    = model.lengths.transexon ;  
  end ;
  if use.segments.intercistronic & use.signals.cleave ,
    A(model.state_ids.cleave(1), model.state_ids.SL1) = model.lengths.intergenictrans ;
    A(model.state_ids.cleave(2), model.state_ids.SL1) = model.lengths.intercistronic ;
    A(model.state_ids.cleave(2), model.state_ids.SL2) = model.lengths.intercistronic ;
  end ;
end ;

% gene start
A(model.state_ids.seq_start,gene_start_states_in) = model.lengths.intergenic ;
%if use.model.use_multi
%  A(model.state_ids.intergenic_long, gene_start_states_in) = model.lengths.intergenic ; 
%end
    
A(gene_start_states_out,model.state_ids.tis(1)) = model.lengths.utr5exon ;
A(gene_start_states_out,model.state_ids.don(1)) = model.lengths.utr5exon ;
% spliced 5'utr
A(model.state_ids.acc(1), model.state_ids.don(1)) = model.lengths.utr5exon ;
A(model.state_ids.don(1), model.state_ids.acc(1)) = model.lengths.intron ;
% last 5'utr exon
A(model.state_ids.acc(1),model.state_ids.tis(1)) = model.lengths.utr5exon ;
A(model.state_ids.don(1), model.state_ids.tis(2)) = model.lengths.intron ; 


%%%%%%%%%%%
% CDS STUFF
%%%%%%%%%%%
A(mRNA_start_states_out, model.state_ids.cdsStop(1)) = model.lengths.single_cds_exon ; % trans(2) ok 
A(mRNA_start_states_out, model.state_ids.don(2:7)) = model.lengths.first_cds_exon ; % trans(2) ok
for i=2:7,
  A(model.state_ids.don(i), model.state_ids.acc(i)) = model.lengths.intron ;
  if isfield(model.state_ids, 'intron_long')
    A(model.state_ids.don(i), model.state_ids.intron_long(i-1)) = model.lengths.intron;
    A(model.state_ids.intron_long(i-1), model.state_ids.acc(i)) = model.lengths.intron_long;
    A(model.state_ids.intron_long(i-1),model.state_ids.intron_long(i-1)) = model.lengths.intron_long;
  end
end ;
A(model.state_ids.acc(2:7), model.state_ids.don(2:7)) = model.lengths.middle_cds_exon ;
A(model.state_ids.acc(2:7), model.state_ids.cdsStop(1)) = model.lengths.last_cds_exon ;
% agstop 
A(model.state_ids.don(2), model.state_ids.cdsStop(2)) = model.lengths.intron ;



%%%%%%%%%%%
% 3' UTR
%%%%%%%%%%%
A(model.state_ids.cdsStop, model.state_ids.don(8)) = model.lengths.utr3exon ;
A(model.state_ids.don(8), model.state_ids.acc(8)) = model.lengths.intron ;
A(model.state_ids.acc(8), model.state_ids.don(8)) = model.lengths.utr3exon ;

if use.signals.polya,
  % polyA has to be in last 3'UTR exon
  A(model.state_ids.cdsStop, model.state_ids.polya) = model.lengths.utr3exon ;
  A(model.state_ids.acc(8), model.state_ids.polya)=model.lengths.utr3exon ;
end ;

gene_stop_states = [];    
if isfield(use.segments, 'rna_seq_polya')&&use.segments.rna_seq_polya
	gene_stop_states = [gene_stop_states model.state_ids.cleave(1)];
	A(model.state_ids.cdsStop, model.state_ids.rna_seq_polya) = model.lengths.utr3exon;
	A(model.state_ids.rna_seq_polya, model.state_ids.cleave) = model.lengths.rna_seq_polya ;
	A(model.state_ids.acc(8), model.state_ids.rna_seq_polya) = model.lengths.utr3exon ;
else
	if use.signals.cleave,
	  gene_stop_states = [gene_stop_states model.state_ids.cleave(1)];
	  A(model.state_ids.cdsStop, model.state_ids.cleave) = model.lengths.utr3exon ;
	  A(model.state_ids.acc(8), model.state_ids.cleave) = model.lengths.utr3exon ;
	  if use.signals.polya,
	    A(model.state_ids.polya, model.state_ids.cleave) = model.lengths.polya_tail ;
	  end
	else
	  gene_stop_states = [gene_stop_states model.state_ids.cdsStop];
	  if use.signals.polya,
	    gene_stop_states = [gene_stop_states model.state_ids.polya];
	  end
	end ;
end

A(gene_stop_states, model.state_ids.seq_end)    = model.lengths.intergenic ; 
%A(gene_stop_states, model.state_ids.intergenic_long) = model.lengths.intergenic ; 
if use.model.use_multi,
  A(gene_stop_states, gene_start_states_in) = model.lengths.intergenic ;    
end ;

if use.segments.intercistronic 
  for k=model.state_ids.SL1,
    idx = find(A(:,k)==model.lengths.intergenic) ;
    A(idx,k) = model.lengths.intergenictrans ;   
  end ;
end
%%%%%%%%%%%
% end transition definition
%%%%%%%%%%%


%%% technical stuff for intergenic_long

transitions = [];
for i=1:length(gene_start_states_in)
  transitions = [transitions; model.state_ids.seq_start gene_start_states_in(i)];
  for j=1:length(gene_stop_states)
    transitions = [transitions; gene_stop_states(j) gene_start_states_in(i)];
    transitions = [transitions; gene_stop_states(j) model.state_ids.seq_end] ;
  end
end

%model.intergenic_long.transitions = unique(transitions,'rows');
%model.intergenic_long.max_short_length = model.lengths_range.intron(2) ;
%model.intergenic_long.step = model.lengths_range.intron(2); % use the intron length as maximal step size
%model.intergenic_long.long_segment_parts = model.lengths_range.intron(2)/2-25 ;
%model.intergenic_long.step = 30 ; % makes about 30 (step) * 40 (grid)=1200nt per step

%%% 

q(gene_stop_states) = 0 ;

model.gene_start_states_in = gene_start_states_in;
model.gene_start_states_out = gene_start_states_out;
model.mRNA_start_states_out = mRNA_start_states_out;
model.gene_stop_states = gene_stop_states;

model.A = A;
model.p = p;
model.q = q;
model.orf_info = orf_info;
model.active_transitions = find(~isinf(A)) ;
model.cnt_transitions = length(model.active_transitions) ;


