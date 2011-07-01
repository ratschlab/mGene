function model = define_states(model)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define states
model.states(1) = struct ;
cnt=0 ;

% sequence start 
cnt=cnt+1 ;
model.states(cnt).consensus={} ;
model.states(cnt).consensus_pos=[] ;
model.states(cnt).non_consensus={} ;
model.states(cnt).non_consensus_pos=[] ;
model.states(cnt).name='seq_start' ;
model.states(cnt).signal=[] ;
model.state_ids.seq_start=cnt ;

signal_names = fieldnames(model.signals) ;
for i=1:model.cnt_signals
  if strcmp(signal_names{i}, 'transacc'), continue; end ;
  if ~strcmp(signal_names{i},'acc') && ~strcmp(signal_names{i},'don') ...
        && ~strcmp(signal_names{i}, 'tis') && ~strcmp(signal_names{i}, 'cdsStop')...
        && ~strcmp(signal_names{i}, 'trans') && ~strcmp(signal_names{i}, 'cleave') 
    cnt=cnt+1 ;
    model.states(cnt).name = signal_names{i};
    model.states(cnt).consensus = getfield(model.consensus,signal_names{i});
    model.states(cnt).consensus_pos = getfield(model.consensus_pos,signal_names{i});
    model.states(cnt).non_consensus = {} ;
    model.states(cnt).non_consensus_pos = [] ;
    model.states(cnt).signal = getfield(model.signals,signal_names{i}) ;
    model.state_ids=setfield(model.state_ids, signal_names{i},cnt) ;
  elseif strcmp(signal_names{i},'cleave') 
    cnt=cnt+1 ;
    model.states(cnt).name = signal_names{i};
    model.states(cnt).consensus = getfield(model.consensus,signal_names{i});
    model.states(cnt).consensus_pos = getfield(model.consensus_pos,signal_names{i});
    model.states(cnt).non_consensus = {} ;
    model.states(cnt).non_consensus_pos = [] ;
    model.states(cnt).signal = getfield(model.signals,signal_names{i}) ;
    model.state_ids=setfield(model.state_ids, signal_names{i},[cnt]) ;
    if model.use.segments.intercistronic
      cnt=cnt+1 ;
      model.states(cnt).name = signal_names{i};
      model.states(cnt).consensus = getfield(model.consensus,signal_names{i});
      model.states(cnt).consensus_pos = getfield(model.consensus_pos,signal_names{i});
      model.states(cnt).non_consensus = {} ;
      model.states(cnt).non_consensus_pos = [] ;
      model.states(cnt).signal = getfield(model.signals,signal_names{i}) ;
      model.state_ids=setfield(model.state_ids, signal_names{i},[cnt-1 cnt]) ;
    end
  elseif strcmp(signal_names{i},'tis') | strcmp(signal_names{i},'cdsStop')
    cnt=cnt+1 ;
    model.states(cnt).name = signal_names{i};
    model.states(cnt).consensus = getfield(model.consensus,signal_names{i});
    model.states(cnt).consensus_pos = getfield(model.consensus_pos,signal_names{i});
    model.states(cnt).non_consensus = {} ;
    model.states(cnt).non_consensus_pos = [] ;
    model.states(cnt).signal = getfield(model.signals,signal_names{i}) ;
    cnt=cnt+1 ;
    model.states(cnt).name = signal_names{i};
    model.states(cnt).non_consensus = {} ;
    model.states(cnt).non_consensus_pos = [] ;
    model.state_ids=setfield(model.state_ids, signal_names{i},[cnt-1 cnt]) ;
    if strcmp(signal_names{i}, 'tis')
    model.states(cnt).name = 'tis_acc' ;
      model.states(cnt).consensus = {'agatg'} ;
      model.states(cnt).consensus_pos = 2 ;
      model.states(cnt).signal = [model.signals.acc model.signals.tis] ;%%%% !!! 2 Zahlen
    else
      model.states(cnt).name = 'cdsStop_acc';
      model.states(cnt).consensus = {'agtaa', 'agtag', 'agtga'} ;
      model.states(cnt).consensus_pos = [2 2 2] ;
      model.states(cnt).signal = [model.signals.acc model.signals.cdsStop] ;%%%% !!! 2 Zahlen
    end
  elseif strcmp(signal_names{i},'trans') & model.use.states.SL1 & model.use.states.SL2
    cnt=cnt+1 ;
    model.states(cnt).name = 'SL1';
    model.states(cnt).consensus = getfield(model.consensus,'trans');
    model.states(cnt).consensus_pos = getfield(model.consensus_pos,'trans');
    model.states(cnt).non_consensus = {} ;
    model.states(cnt).non_consensus_pos = [] ;
    model.states(cnt).signal = [model.signals.trans model.signals.transacc]  ;
    cnt=cnt+1 ;
    model.states(cnt).name = 'SL1';
    model.states(cnt).consensus = {'agatg'} ; 
    model.states(cnt).consensus_pos = 2;
    model.states(cnt).non_consensus = {} ;
    model.states(cnt).non_consensus_pos = [] ;
    model.states(cnt).signal = [model.signals.trans model.signals.transacc model.signals.tis]  ;
    model.state_ids = setfield(model.state_ids, 'SL1', [cnt-1 cnt]) ;
    cnt=cnt+1 ;
    model.states(cnt).name = 'SL2';
    model.states(cnt).consensus = getfield(model.consensus,'trans');
    model.states(cnt).consensus_pos = getfield(model.consensus_pos,'trans');
    model.states(cnt).non_consensus = {} ;
    model.states(cnt).non_consensus_pos = [] ;
    model.states(cnt).signal = [model.signals.trans model.signals.transacc]  ;
    cnt=cnt+1 ;
    model.states(cnt).name = 'SL2';
    model.states(cnt).consensus = {'agatg'} ; 
    model.states(cnt).consensus_pos = 2;
    model.states(cnt).non_consensus = {} ;
    model.states(cnt).non_consensus_pos = [] ;
    model.states(cnt).signal = [model.signals.trans model.signals.transacc model.signals.tis]  ;%%%% !!! 3 Zahlen
    model.state_ids = setfield(model.state_ids, 'SL2', [cnt-1 cnt]) ;
  elseif strcmp(signal_names{i},'trans') & ~model.use.states.SL1 & ~model.use.states.SL2
    cnt=cnt+1 ;
    model.states(cnt).name = 'trans';
    model.states(cnt).consensus = getfield(model.consensus,'trans');
    model.states(cnt).consensus_pos = getfield(model.consensus_pos,'trans');
    model.states(cnt).non_consensus = {} ;
    model.states(cnt).non_consensus_pos = [] ;
    model.states(cnt).signal = [model.signals.trans model.signals.transacc]  ;
    cnt=cnt+1 ;
    model.states(cnt).name = 'trans';
    model.states(cnt).consensus = {'agatg'} ; 
    model.states(cnt).consensus_pos = 2;
    model.states(cnt).non_consensus = {} ;
    model.states(cnt).non_consensus_pos = [] ;
    model.states(cnt).signal = [model.signals.trans model.signals.transacc model.signals.tis]  ;
    model.state_ids = setfield(model.state_ids, 'SL1', [cnt-1 cnt]) ;
  elseif strcmp(signal_names{i},'don')
    cnt_don = cnt ;
    cnt=cnt+8 ;
    for id = cnt_don+1:cnt_don+8 ;
      model.states(id).consensus = getfield(model.consensus,signal_names{i});
      model.states(id).consensus_pos = getfield(model.consensus_pos,signal_names{i});
      model.states(id).non_consensus = {} ;
      model.states(id).non_consensus_pos = [] ;
      model.states(id).signal=model.signals.don ;
    end ;
    model.states(cnt_don+1).name='UTR5_don' ;
    model.states(cnt_don+2).name='CDS_don_0' ;
    model.states(cnt_don+3).name='CDS_don_1a' ;
    model.states(cnt_don+3).non_consensus = {'t'} ;
    model.states(cnt_don+3).non_consensus_pos= [1] ;
    model.states(cnt_don+4).name='CDS_don_1b' ;
    model.states(cnt_don+4).consensus = {'tgt','tgc'} ;
    model.states(cnt_don+4).consensus_pos= [1,1] ;
    model.states(cnt_don+5).name='CDS_don_2a' ;
    model.states(cnt_don+5).non_consensus = {'ta', 'tg'} ;
    model.states(cnt_don+5).non_consensus_pos= [2,2] ;
    model.states(cnt_don+6).name='CDS_don_2b' ;
    model.states(cnt_don+6).consensus = {'tagt', 'tagc'} ;
    model.states(cnt_don+6).consensus_pos= [2,2] ;
    model.states(cnt_don+7).name='CDS_don_2c' ;
    model.states(cnt_don+7).consensus = {'tggt', 'tggc'} ;
    model.states(cnt_don+7).consensus_pos= [2,2] ;
    model.states(cnt_don+8).name='UTR3_don' ;
    
    model.state_ids.don=cnt_don+1:cnt_don+8 ;
  elseif strcmp(signal_names{i},'acc')
    cnt_acc = cnt ;
    cnt=cnt+8 ;
    for id = cnt_acc+1:cnt_acc+8 ;
      model.states(id).consensus = {'ag'} ;
      model.states(id).consensus_pos= [2] ;
      model.states(id).non_consensus={} ;
      model.states(id).non_consensus_pos=[] ;
      model.states(id).signal=model.signals.acc ;
    end ;
    model.states(cnt_acc+1).name='UTR5_acc' ;
    model.states(cnt_acc+2).name='CDS_acc_0' ;
    model.states(cnt_acc+3).name='CDS_acc_1a' ;
    model.states(cnt_acc+4).name='CDS_acc_1b' ;
    model.states(cnt_acc+4).non_consensus = {'aa','ag','ga'} ;
    model.states(cnt_acc+4).non_consensus_pos= [0,0,0] ;
    model.states(cnt_acc+5).name='CDS_acc_2a' ;
    model.states(cnt_acc+6).name='CDS_acc_2b' ;
    model.states(cnt_acc+6).consensus = {'agt', 'agc'} ;
    model.states(cnt_acc+6).consensus_pos= [2, 2] ;
    model.states(cnt_acc+7).name='CDS_acc_2c' ;
    model.states(cnt_acc+7).non_consensus = {'a'} ;
    model.states(cnt_acc+7).non_consensus_pos= [0] ;
    model.states(cnt_acc+8).name='UTR3_acc' ;
 
    model.state_ids.acc=cnt_acc+1:cnt_acc+8 ;

    if isfield(model.use, 'intron_long_state')&& model.use.intron_long_state
      cnt_intron_long = cnt;
      cnt = cnt+6;
      model.states(cnt_intron_long+1).consensus={} ;
      model.states(cnt_intron_long+1).consensus_pos=[] ;
      model.states(cnt_intron_long+1).non_consensus={} ;
      model.states(cnt_intron_long+1).non_consensus_pos=[] ;
      model.states(cnt_intron_long+1).name='intron_long_0' ;
      model.states(cnt_intron_long+1).signal=[] ;

      model.states(cnt_intron_long+2).consensus={} ;
      model.states(cnt_intron_long+2).consensus_pos=[] ;
      model.states(cnt_intron_long+2).non_consensus={} ;
      model.states(cnt_intron_long+2).non_consensus_pos=[] ;
      model.states(cnt_intron_long+2).name='intron_long_1a' ;
      model.states(cnt_intron_long+2).signal=[] ;
    
      model.states(cnt_intron_long+3).consensus={} ;
      model.states(cnt_intron_long+3).consensus_pos=[] ;
      model.states(cnt_intron_long+3).non_consensus={} ;
      model.states(cnt_intron_long+3).non_consensus_pos=[] ;
      model.states(cnt_intron_long+3).name='intron_long_1b' ;
      model.states(cnt_intron_long+3).signal=[] ;

      model.states(cnt_intron_long+4).consensus={} ;
      model.states(cnt_intron_long+4).consensus_pos=[] ;
      model.states(cnt_intron_long+4).non_consensus={} ;
      model.states(cnt_intron_long+4).non_consensus_pos=[] ;
      model.states(cnt_intron_long+4).name='intron_long_2a' ;
      model.states(cnt_intron_long+4).signal=[] ;

      model.states(cnt_intron_long+5).consensus={} ;
      model.states(cnt_intron_long+5).consensus_pos=[] ;
      model.states(cnt_intron_long+5).non_consensus={} ;
      model.states(cnt_intron_long+5).non_consensus_pos=[] ;
      model.states(cnt_intron_long+5).name='intron_long_2b' ;
      model.states(cnt_intron_long+5).signal=[] ;

      model.states(cnt_intron_long+6).consensus={} ;
      model.states(cnt_intron_long+6).consensus_pos=[] ;
      model.states(cnt_intron_long+6).non_consensus={} ;
      model.states(cnt_intron_long+6).non_consensus_pos=[] ;
      model.states(cnt_intron_long+6).name='intron_long_1a' ;
      model.states(cnt_intron_long+6).signal=[] ;

      model.state_ids.intergenic_long= cnt_intron_long+1:cnt_intron_long+6;
    end
%    cnt
  end
end

% intergenic_long
%cnt=cnt+1 ;
%model.states(cnt).consensus={} ;
%model.states(cnt).consensus_pos=[] ;
%model.states(cnt).non_consensus={} ;
%model.states(cnt).non_consensus_pos=[] ;
%model.states(cnt).name='intergenic_long' ;
%model.states(cnt).signal=[] ;
%model.state_ids.intergenic_long=cnt ;

if isfield(model.use.contents, 'rna_seq_polya')&&model.use.contents.rna_seq_polya
	cnt=cnt+1 ;
	model.states(cnt).consensus={} ;
	model.states(cnt).consensus_pos=[] ;
	model.states(cnt).non_consensus={} ;
	model.states(cnt).non_consensus_pos=[] ;
	model.states(cnt).name='rna_seq_polya' ;
	model.states(cnt).signal=[] ;
	model.state_ids.rna_seq_polya = cnt ;
end

if isfield(model.use, 'non_coding')&&model.use.non_coding
	cnt=cnt+1 ;
	model.states(cnt).consensus={} ;
	model.states(cnt).consensus_pos=[] ;
	model.states(cnt).non_consensus={} ;
	model.states(cnt).non_consensus_pos=[] ;
	model.states(cnt).name='tss_nc' ;
	model.states(cnt).signal= getfield(model.signals, 'tss') ;
	model.state_ids.tss_nc = cnt ;
	cnt=cnt+1 ;
    model.states(cnt).consensus = getfield(model.consensus, 'acc');
    model.states(cnt).consensus_pos = getfield(model.consensus_pos, 'acc');
	model.states(cnt).non_consensus={} ;
	model.states(cnt).non_consensus_pos=[] ;
	model.states(cnt).name='acc_nc' ;
	model.states(cnt).signal= getfield(model.signals, 'acc') ;
	model.state_ids.acc_nc = cnt ;
	cnt=cnt+1 ;
    model.states(cnt).consensus = getfield(model.consensus, 'don');
    model.states(cnt).consensus_pos = getfield(model.consensus_pos, 'don');
	model.states(cnt).non_consensus={} ;
	model.states(cnt).non_consensus_pos=[] ;
	model.states(cnt).name='don_nc' ;
	model.states(cnt).signal= getfield(model.signals, 'don') ;
	model.state_ids.don_nc = cnt ;
	cnt=cnt+1 ;
	model.states(cnt).consensus={} ;
	model.states(cnt).consensus_pos=[] ;
	model.states(cnt).non_consensus={} ;
	model.states(cnt).non_consensus_pos=[] ;
	model.states(cnt).name='cleave_nc' ;
	model.states(cnt).signal= getfield(model.signals, 'cleave') ;
	model.state_ids.cleave_nc = cnt ;
end

% sequence end
cnt=cnt+1 ;
model.states(cnt).consensus={} ;
model.states(cnt).consensus_pos=[] ;
model.states(cnt).non_consensus={} ;
model.states(cnt).non_consensus_pos=[] ;
model.states(cnt).name='seq_end' ;
model.states(cnt).signal=[] ;
model.state_ids.seq_end = cnt ;

model.cnt_states = cnt ; 
