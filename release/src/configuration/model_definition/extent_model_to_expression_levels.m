function extended_model = extent_model_to_expression_levels(model, num_levels)

if nargin==0
	%load /fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/human/lsl/abinitio_nuc_pos/output/lsl/data/training_PAR.mat
	load /fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/elegans/lsl/rna_seq_polya/output/lsl/data/training_PAR.mat
	model = PAR.model
	num_levels = 3;
end


% allocate space for the new transition matrices
% and insert the old matrix along the main diagonal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tp_dim3 = size(model.transition_pointers, 3);
num_states = size(model.A, 1);
A = -inf(num_levels*size(model.A));
tp = zeros(num_levels*num_states, num_levels*num_states, tp_dim3);
for j = 0:num_levels-1
	last_idx = j*num_states+1;
	current_idx = (j+1)*num_states;
	A(last_idx:current_idx, last_idx:current_idx) = model.A;
	for k = 1:tp_dim3
		tp(last_idx:current_idx, last_idx:current_idx, k) = model.transition_pointers(:,:,k)';
	end
end


% allow transitons in intergenic regions between all 
% expression levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp = sparse(model.transition_pointers(:,:,2)'==model.contents.intergenic)

[rows, cols, tmp] = find(sp);

for j = 1:num_levels
	for k = setdiff(1:num_levels, j)
		for r = 1:length(rows)
			row = (j-1)*num_states+rows(r);
			col = (k-1)*num_states+cols(r);
			A(row, col) = 10;
			%A(col, row) = 10;
			for l = 1:tp_dim3
				tp(row, col,l) = model.transition_pointers(rows(r), cols(r),l);
			end
		end
	end
end

%% duplicate states and orf info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
states = [];
orf_info = [];
for j = 1:num_levels
	new_states = model.states;
	for k = 1:num_states
		new_states(k).name = sprintf('%s_lev%i',model.states(k).name, j);
	end
	states = [states new_states];
	orf_info = [orf_info; model.orf_info];
end


%% remove seq_start and seq_end states between expression levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
keep_idx = [];
for j = 2:length(states)-1%%first and last seq_start/seq_end state are needed
	if isempty(strfind(states(j).name, 'seq_end'))&&isempty(strfind(states(j).name, 'seq_start'))
		keep_idx = [keep_idx j];
	end
end
keep_idx = [1 keep_idx num_levels*num_states];
A = A(keep_idx, keep_idx);
states = states(keep_idx);
transition_pointers = zeros(length(keep_idx), length(keep_idx), tp_dim3);
for l = 1:tp_dim3
	transition_pointers(:,:,l) = tp(keep_idx,keep_idx,l)';% transpose
end
orf_info = orf_info(keep_idx, :);

% rename seq_start and seq_end again
states(1).name = model.states(1).name;
states(end).name = model.states(end).name;


% map states to the signals
old_states = fieldnames(model.state_ids);
state_ids = struct;
for j = 1:length(states)
	name = states(j).name;
	field = '';
	for k = 1:length(old_states)
		if strfind(name, old_states{k})
			field = old_states{k};
			break;
		end
	end
	if isempty(field)
		error('could not map state\n');
	end
	if isfield(state_ids, field)
		state_ids.(field) = [state_ids.(field), j];
	else
		state_ids.(field) = j;
	end
end

% start and end state distribution
p = -inf(length(states), 1);
p(state_ids.seq_start) = 0;
q = -inf(length(states), 1);
q(state_ids.seq_end) = 0;

a_trans = define_a_trans(A, transition_pointers, model.seg_links);

% assign to model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extended_model = model;
extended_model.states = states;
extended_model.cnt_states = length(states);
extended_model.orf_info = orf_info;
extended_model.state_ids = state_ids;
extended_model.A = A;
extended_model.active_transitions = find(~isinf(A));
extended_model.cnt_transitions = length(extended_model.active_transitions);
extended_model.transition_pointers = transition_pointers;
extended_model.a_trans = a_trans;
extended_model.p = p;
extended_model.q = q;
extended_model.cnt_parameters = extended_model.cnt_transitions + model.cnt_plifs*model.bins ;
extended_model = plif_to_param_mapping(extended_model);
extended_model.use.expression_levels = num_levels;

print_graph_dot(extended_model, '~/tmp/model_extended', 0)



