function gene = apply_transcript_idx(gene, idx, fields, do_update_start_stop);
% gene = apply_transcript_idx(gene, idx, fields, do_update_start_stop);

if nargin<4,
    do_update_start_stop=0 ;
end ;

if nargin<3 || isequal(fields, []),
    fields = default_gene_field_names() ;
end ;

for i=1:length(fields)
	if isfield(gene, fields{i}) && (isempty(idx) || length(gene.(fields{i}))>=max(idx))
		gene.(fields{i}) 	= gene.(fields{i})(idx);
	end
end ;

if do_update_start_stop,
    gene = update_start_stop(gene);
end ;
