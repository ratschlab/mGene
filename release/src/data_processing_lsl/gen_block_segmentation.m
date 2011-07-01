function blocks = gen_block_segmentation(genes,blocks,model)
% blocks = gen_block_segmentation(genes,blocks,model)
%  
% genes: need to have a field 'segments' with a gene segmentation in global chr-coordinates (allways on positive strand)
% 
% segmentation: 4 columns
%         col 1: start of segment: must be equal to col2 of previous line
%         col 2: end of segment
%         col 3: segment type (e.g. intergenic, exon, intron,...) integer coded. definition: model.segments 
%         col 4: gene count: 1:length(blocks.gene_ids); 0 if segment is intergenic

use_weights = isfield(genes, 'transcript_weights');

for b=1:length(blocks)
  current = 1;
  seg = {[1 length(blocks(b).seq) model.segments.intergenic 0]};
  if use_weights
	weights = -inf;
	assert(length(blocks(b).gene_ids)==1)
  end
  for g=1:length(blocks(b).gene_ids)
    gene = genes(blocks(b).gene_ids(g));
    if ~gene.is_complete==1, continue; end
    if gene.start<blocks(b).start||gene.start>blocks(b).stop||gene.stop<blocks(b).start||gene.stop>blocks(b).stop
      error('gene boundaries not within block boundaries')
    end
    assert(gene.id==blocks(b).gene_ids(g))
    for s = 1:length(seg)
      for t=1:length(gene.exons)
		%% filter transcripts
        if length(gene.transcript_valid)>=t&&~gene.transcript_valid(t), 
			continue; 
		end
		if length(gene.transcript_complete)>=t&&~gene.transcript_complete(t)
			continue;
		end
		if isfield(gene, 'transcript_coding') && gene.transcript_coding(t)~=1
			continue
		end

        if blocks(b).strand=='+'%
          local_segments = [gene.segments{t}(:,1:2)-blocks(b).start+1 gene.segments{t}(:,3)];
          local_segments(:,end+1)=g;
        else %
          local_segments = [blocks(b).stop-gene.segments{t}(end:-1:1,2:-1:1)+1 gene.segments{t}(end:-1:1,3)];
          local_segments(:,end+1)=length(blocks(b).gene_ids)-g+1;
        end%
        seg{end+1} = insert_gene_segments(seg{s},local_segments,model,blocks(b).strand);
		if use_weights
			weights(end+1) = gene.transcript_weights(t);
		end
      end
    end
    seg(1:s)=[];
	if use_weights
		weights(1:s)= [];
	end
  end
  if use_weights
    assert(length(weights)==length(seg));
  end
  for j=1:length(seg)
    blocks(b).truth(j).segments=seg{j};
	if use_weights
      blocks(b).truth(j).transcript_weights=weights(j);
	end
  end
end
return

function segments = insert_gene_segments(segments,gene_segments,model,strand)
  % example for positive strand: 
  % segments = [1 1000 0]
  % gene_segments = [300 400 1]
  %                 [400 500 2]
  % result = [1   300  0]
  %          [300 400  1]
  %          [400 500  2]
  %          [500 1000 0]
  %
if gene_segments(1,1)==1
	% nasty special case
	gene_segments(1,1)=2;
end
line = find(segments(:,1)<gene_segments(1,1),1,'last');

% gene_segments have to fit into one single intergenic region
if ~(segments(line,3)==model.segments.intergenic)
  segments = [];
  return
end
assert(segments(line,2)>=gene_segments(end,2))

new_line = [gene_segments(end,2) segments(line,2)  model.segments.intergenic 0];
segments(line,:)=[segments(line,1) gene_segments(1,1) model.segments.intergenic 0];
segments = [segments(1:line,:); gene_segments; new_line; segments(line+1:end,:)]; 
return


