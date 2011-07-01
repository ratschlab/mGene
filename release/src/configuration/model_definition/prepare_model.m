function [model,lengths] = prepare_model(PAR, blocks)
% [model,lengths] = prepare_model(PAR, blocks)
% - getting boundaries for plifs
  
model = PAR.model; 
%fn_contents = PAR.FN.input_lsl.fn_content
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boundaries for signal and content detectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0%fexist(PAR.FN.input_lsl.fn_boundary_model),
  fprintf('using boundary model in %s\n', PAR.FN.input_lsl.fn_boundary_model) ;
  L = load(PAR.FN.input_lsl.fn_boundary_model) ;
  model.boundaries = L.boundaries ;  
  return
end

if nargin<2
  blocks = load_struct(PAR.FN.input_lsl.fn_training_blocks, 'blocks') ;  
end
INF = 1e20 ; 


%%% signal detectors
  
boundaries = struct;
boundaries.signals = struct;
boundaries.contents = struct;
boundaries.lengths = struct;

sig_names = fieldnames(model.signals) ;
for s = 1:length(sig_names)
  range = getfield(model.signals_range, sig_names{s}) ;
  boundaries.signals = setfield(boundaries.signals, sig_names{s},[-INF linspace(range(1),range(2),model.bins-1) INF]);
end

%%% content detectors

content_names = fieldnames(model.contents) ;
for s = 1:length(content_names)
  range = getfield(model.contents_range, content_names{s}) ;
  boundaries.contents = setfield(boundaries.contents, content_names{s},[-INF linspace(range(1),range(2),model.bins-1) INF]);
end

for j = 1:length(model.track_names)
  fieldname = sprintf('track_%i',j);
  boundaries.(fieldname) = struct;
  names = fieldnames(model.(fieldname)) ;
  range = determine_range(j, blocks);
  for s = 1:length(names)
    boundaries.(fieldname) = setfield(boundaries.(fieldname), names{s},[-INF linspace(range(1),range(2),model.bins-1) INF]);
  end
end
for j = 1:length(model.segment_feature_names);
  fieldname = sprintf('segment_feature_%i',j);
  boundaries.(fieldname) = struct;
  names = fieldnames(model.(fieldname)) ;
  range = [0 18] ; 
  for s = 1:length(names)
    boundaries.(fieldname) = setfield(boundaries.(fieldname), names{s},[-INF linspace(range(1),range(2),model.bins-1) INF]);
  end
end
for j = 1:length(model.segment_feature_names);
  fieldname = sprintf('segment_score_%i',j);
  boundaries.(fieldname) = struct;
  names = fieldnames(model.(fieldname)) ;
  range = [1 30]; 
  for s = 1:length(names)
    boundaries.(fieldname) = setfield(boundaries.(fieldname), names{s},[-INF linspace(range(1),range(2),model.bins-1) INF]);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% length distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

length_names = fieldnames(model.lengths) ;
lengths = [] ;
for s = 1:model.cnt_lengths
  lengths = setfield(lengths, length_names{s}, []) ;
end

segments_names = fieldnames(model.segments) ;
for id = 1:length(blocks),
  block = blocks(id) ;
  segments = block.truth(1).segments ;
  idx_genes=setdiff(unique(segments(:,4)),0);
  for i=1:length(idx_genes)
    genes{i} = segments(segments(:,4)==idx_genes(i),1:3);
  end
  for s = 1:length(segments_names)
    s_id = getfield(model.segments, segments_names{s}) ;
    if s_id == model.segments.cds_exon,
      for j=1:length(genes),
        gene_segments = genes{j} ;
        idx = find(gene_segments(:,3)==s_id) ;
        if length(idx)==1,
          lengths.single_cds_exon = [lengths.single_cds_exon ;gene_segments(idx,2)-gene_segments(idx,1)];
        elseif length(idx)>1,
          lengths.first_cds_exon = [lengths.first_cds_exon ;gene_segments(idx(1),2)-gene_segments(idx(1),1)] ;
          lengths.last_cds_exon = [lengths.last_cds_exon ;gene_segments(idx(end),2)-gene_segments(idx(end),1)] ;
          lengths.middle_cds_exon = [lengths.middle_cds_exon ;gene_segments(idx(2:end-1),2)-gene_segments(idx(2:end-1),1)] ;
        else
          % this is a non coding gene
          %assert(0) ;
        end ;
      end ;
	elseif isfield(model.segments, 'nc_exon') && s_id == model.segments.nc_exon
		gene_segments = genes{j};
		idx = find(gene_segments(:,3)==s_id);
		if length(idx)==1
			lengths.middle_nc_exon = [lengths.middle_nc_exon; gene_segments(idx,2)-gene_segments(idx,1)] ;
		elseif length(idx)==2
			lengths.first_nc_exon = [lengths.first_nc_exon; gene_segments(idx(1),2)-gene_segments(idx(1),1)] ;
			lengths.last_nc_exon = [lengths.last_nc_exon; gene_segments(idx(2),2)-gene_segments(idx(2),1)] ;
		elseif length(idx)>2
			lengths.first_nc_exon = [lengths.first_nc_exon; gene_segments(idx(1),2)-gene_segments(idx(1),1)] ;
			lengths.middle_nc_exon = [lengths.middle_nc_exon; gene_segments(idx(2:end-1),2)-gene_segments(idx(2:end-1),1)] ;
			lengths.last_nc_exon = [lengths.last_nc_exon; gene_segments(idx(end),2)-gene_segments(idx(end),1)] ;
		end
    else
      idx = find(segments(:,3)==s_id) ;
      list = getfield(lengths, segments_names{s}) ;
      list = [list; segments(idx,2)-segments(idx,1)] ;
      lengths = setfield(lengths, segments_names{s}, list) ;
    end ;
  end
end
for s = 1:length(length_names)
  if isempty(lengths.(length_names{s}))
    warning('nonexisting field %s. Using default plif values',length_names{s});
    sout = log(logspace(log10(5),log10(1000),1000)+1);
    boundaries.lengths = setfield(boundaries.lengths,length_names{s}, [-INF sout(ceil([1/model.bins:1/model.bins:1-1/model.bins]*length(sout))) INF]) ;
    continue;
  end
  %range = getfield(model.lengths_range, length_names{s}) ;
  sout = getfield(lengths,length_names{s}) ;
  %sout(sout<range(1)|sout>range(2)) = [];
  %sout = sort(log(sout+1))' ;
  %fprintf(1, '%s max_val: %i\n', length_names{s}, max(sout))
  sout = log(sout+1)' ;
  %if length(sout)<100, 
    %fprintf('too few points: %s %i\n',length_names{s}, length(sout));  
    %sout = log(logspace(log10(5),log10(1000),1000)+1) ; 
    %sout = log(logspace(log10(range(1)),log10(range(2)),1000)+1) ; 
  %end ;
  %boundaries.lengths = setfield(boundaries.lengths,length_names{s}, [-INF sout(ceil([1/model.bins:1/model.bins:1-1/model.bins]*length(sout))) INF]) ;
  %keyboard
  boundaries.lengths = setfield(boundaries.lengths,length_names{s}, [-INF prctile(sout,(1/(model.bins-1):1/(model.bins-1):1)*100) INF]) ;
end
model.boundaries=boundaries ;


fprintf('saving boundary model in %s\n', PAR.FN.input_lsl.fn_boundary_model) ;
boundaries  = model.boundaries;
save(PAR.FN.input_lsl.fn_boundary_model, 'boundaries', '-v7') ;

return

function range = determine_range(idx, blocks)
  vec = [];
  for j = 1:length(blocks)
    local_pos = find(~isnan(blocks(j).tracks(idx,:)));
    vec = [vec blocks(j).tracks(idx, local_pos)];
  end
  range = [min(vec), max(vec)];
  %range = my_prctile(vec, [1 99]);
return
