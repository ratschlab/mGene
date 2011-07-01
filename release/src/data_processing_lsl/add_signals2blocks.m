function blocks = add_signals2blocks(blocks, PAR, signal_names, which_dir, offset,which_fields)
% blocks = add_signals2blocks(blocks, PAR, signal_names, which_dir, offset)

if nargin<5
  warning('I think this should be a parameter')
  offset = 100;  
end
if nargin<4 || isempty(which_dir)
  which_dir = 'pred_dir' ;
end
genome_info = init_genome(PAR.FN.input.fn_genome_config);
which_pos.local = 1;
which_pos.contig = 1;
base_dirs    = {};

for i=1:length(signal_names),
  if isequal(which_dir,'pred_dir')'
    tmp = PAR.FN.output_sig.(signal_names{i});
    base_dirs{i} = tmp.fn_pred;
    fields = {'Conf_cum'};
  elseif isequal(which_dir,'cand_dir')'
    tmp = PAR.FN.input_sig.(signal_names{i});
    base_dirs{i} = tmp.fn_candsites;
    fields = {'label','altgenic','alt_intron','est_conf','est_conf_skip','cDNA_conf','cDNA_conf_skip'};
  end
end
if nargin>5
  fields = which_fields;
end

if ~isfield(blocks,'Signals')
  blocks(1).Signals=struct;
end

% blocks = load_sequence(blocks);
fprintf('Processing %i blocks: ', length(blocks));
for r=1:length(blocks)
  block=blocks(r);
  if mod(r,50)==0, fprintf('%i .. ', r); end ;
  % fill the field genestr
  for s=1:length(signal_names),
    f_name = sprintf('%scontig_%i%s.pos', base_dirs{s},block.chr_num,block.strand);
    d = dir(f_name) ;
	if isempty(base_dirs{s})
		% use zeros instead of signal prediction
		d = [];
		d.bytes = 1;
	end
    if isempty(d), error('file not found: %s',f_name);end
    if d.bytes>0
      block = retrieve_signals(block, base_dirs{s},signal_names{s}, fields, which_pos, offset);
      % block = retrieve_signals(block, base_dirs{s},signal_names{s}, fields,which_pos,0);
    else
      warning('hack!\n')
      block.Signals.(signal_names{s})=struct ;
      block.Signals.(signal_names{s}).pos = [];
      block.Signals.(signal_names{s}).contig_pos = [];
      for k=1:length(fields)
        block.Signals.(signal_names{s}).(fields{k}) =  [];
      end ;
    end
  end
  blocks(r)=block;
end
fprintf('Done.\n') ;

% load /fml/ag-raetsch/share/projects/nGASP/LSL/test_data/cat_1/blocks_2006-12-19.mat

% eof
