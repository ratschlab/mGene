function model = define_signals(model,polya_sixmers)
  
% declare signal types
cnt=0;
signal_start = model.cnt_plifs+1;

if model.use.signals.tss
  model.signals.tss = signal_start+cnt;
  model.consensus.tss = {} ;
  model.consensus_pos.tss = [];
  model.signals_range.tss = [0,1] ;
  cnt=cnt+1;
end

model.signals.tis = signal_start+cnt;
model.consensus.tis = {'atg'} ;
model.consensus_pos.tis = [0];
model.signals_range.tis = [0,1] ;
cnt=cnt+1;

model.signals.don = signal_start+cnt;
model.consensus.don = {'gt','gc'} ;
model.consensus_pos.don = [0,0];
model.signals_range.don = [0,1] ;
cnt=cnt+1;

model.signals.acc = signal_start+cnt;
model.consensus.acc = {'ag'} ;
model.consensus_pos.acc = [2];
model.signals_range.acc = [0,1] ;
cnt=cnt+1;

model.signals.cdsStop = signal_start+cnt ;
model.consensus.cdsStop = {'taa','tga','tag'} ;
model.consensus_pos.cdsStop = [0,0,0];
model.signals_range.cdsStop = [0,1] ;
cnt=cnt+1;

if model.use.signals.polya
  model.signals.polya = signal_start+cnt ;
  model.consensus.polya = polya_sixmers;
  model.consensus_pos.polya = zeros(1,length(polya_sixmers)) ;
  model.signals_range.polya = [0,1] ;
  cnt=cnt+1;
end
if model.use.signals.cleave
  model.signals.cleave = signal_start+cnt ;
  model.consensus.cleave = {} ;
  model.consensus_pos.cleave = [];
  model.signals_range.cleave = [0,1] ; 
  cnt=cnt+1;
end

if model.use.signals.trans
  model.signals.trans = signal_start+cnt;
  model.consensus.trans = {'ag'} ;
  model.consensus_pos.trans = [2];
  model.signals_range.trans = [0,1] ;
  cnt=cnt+1;
end
if model.use.signals.transacc
  model.signals.transacc = signal_start+cnt;
  model.consensus.transacc = {'ag'} ;
  model.consensus_pos.transacc = [2];
  model.signals_range.transacc = [0,1] ;
  cnt=cnt+1;
end
%if isfield(model.use.signals, 'rna_seq_polya')&&model.use.signals.rna_seq_polya
%	model.signals.rna_seq_polya = signal_start+cnt;
%	model.consensus.rna_seq_polya = {} ;
%	model.consensus_pos.rna_seq_polya = [];
%	model.signals_range.rna_seq_polya = [0,1] ;
%	cnt=cnt+1;
%end
model.cnt_signals = length(fieldnames(model.signals));

model.signals_fields = fieldnames(model.signals) ;

model.cnt_plifs = model.cnt_plifs + model.cnt_signals;
