function [SEG, LABEL] = get_cand_contents(P)
% [SEG, LABEL] = get_cand_contents(P)


genes = P.genes ;
dist_to_genes = P.intergenic_dist_to_genes;
use_frame = P.use_frame;
which_segments = P.which_segments;
segment_ids = P.segment_ids;
genome_info = P.genome_info;

clear P

%i=min(find([genes.transcript_complete]));
pos_segments_all = ones(1e6, 1) ;%size(genes(i).segments{1},2));
neg_segments_all = ones(1e6, 1) ;%size(genes(i).segments{1},2));
%i=min(find([genes.transcript_complete]));
%pos_segments_all = ones(1e6,size(genes(i).segments{1},2));
%neg_segments_all = ones(1e6,size(genes(i).segments{1},2));
num_pos = 0;
num_neg = 0;

last_gene_stop=1;
[tmp,idx] = sort([genes.start]);
genes = genes(idx);
for i=1:length(genes)
  gene = genes(i);
  if ~gene.is_valid || ~any(gene.transcript_complete)
    last_gene_stop = gene.stop;
    continue
  end
  % find common/consens segments from different truths
  %----------------------------------------------------------------
  if isempty(gene.segments)
    last_gene_stop = gene.stop;
    continue
  end
  if ~use_frame
    cons_segments = gene.segments{1};
    for t = 2:length(gene.segments)
      if isempty(cons_segments) 
        cons_segments=gene.segments{t} ;
        %break
      end
      if ~gene.transcript_complete(t),
        continue ;
      end ;
      if isempty(gene.segments{t})
        continue
      end
      [tmp,idx1,idx2] = intersect(cons_segments(:,1),gene.segments{t}(:,1) ) ;
      ii = find(cons_segments(idx1,2)==gene.segments{t}(idx2,2));
      ii2 = find(cons_segments(idx1(ii),3)==gene.segments{t}(idx2(ii),3));
      cons_segments = cons_segments(idx1(ii(ii2)),:);
      % doesn't work in octave
      % cons_segments = intersect(cons_segments,gene.segments{t},'rows');
    end
    if isempty(cons_segments)
      cons_segments=zeros(0,3);
    end
    
    % extract lines from common/consens segments 
    %----------------------------------------------------------------
    % if any(cons_segments(2:end,2)-cons_segments(1:end-1,2)<0)
    %   keyboard
    % end
    
    pos_segments = cons_segments(find(ismember(cons_segments(:,3), which_segments{1})), :);
    neg_segments = cons_segments(find(ismember(cons_segments(:,3), which_segments{2})), :);
    
    if ismember(segment_ids.intergenic,which_segments{1})
      assert(isempty(pos_segments))
      if gene.start - dist_to_genes(1)>last_gene_stop
        pos_segments = [pos_segments;  max(last_gene_stop,gene.start-dist_to_genes(2)), gene.start segment_ids.intergenic];
      end
    end

    if ismember(segment_ids.intergenic,which_segments{2})
      if gene.start - dist_to_genes(1)>last_gene_stop
        neg_segments = [neg_segments; max(last_gene_stop,gene.start-dist_to_genes(2)), gene.start segment_ids.intergenic];
      end
    end
  else % use_frame
    cds_exons = gene.cds_exons{1};
    % just checks
    if isempty(cds_exons), continue; end
    if genes(i).strand=='+'
      cds_exons(:,2) = cds_exons(:,2)-1;
    else
      cds_exons(:,1) = cds_exons(:,1)+1;
    end
    [cds,valid_splice] =  load_genomic(gene.chr, char(gene.strand),cds_exons(:,1),cds_exons(:,2),genome_info,1);
    f = find_open_frames(cds,1) ;
    too_short = cds_exons(1,2)-cds_exons(1,1)+1<3||cds_exons(end,2)-cds_exons(end,1)+1<3;
    if mod(length(cds),3)~=0 || f(1)~=1 || isempty(cds) || any(cds(1:3)~='atg') || ~any(ismember({'taa', 'tag', 'tga'}, cds(end-2:end)))|| too_short
      warning('invalid transcript') ;
      continue
    end ;
    % create positive examples
    cds_exons = gene.cds_exons{1};
    %% remove ATG and TAG
    cds_exons(1,1) = cds_exons(1,1)+3;
    cds_exons(end,2) = cds_exons(end,2)-3;
    if genes(i).strand=='+'
      cds_exons(cds_exons(:,1)==1,1) = cds_exons(cds_exons(:,1)==1,1)+2;
      cds_exons(cds_exons(:,1)==1,2) = cds_exons(cds_exons(:,2)==2,1)+1;
    else
      cds_exons(cds_exons(:,1)==1,2) = cds_exons(cds_exons(:,1)==1,2)-2;
      cds_exons(cds_exons(:,2)==2,2) = cds_exons(cds_exons(:,2)==2,2)-1;
    end
    pos_segments = cds_exons(:,1:2);
    %% create negative examples
    cds_exons = gene.cds_exons{1};
    %% remove ATG and TAG
    cds_exons(1,1) = cds_exons(1,1)+3;
    cds_exons(end,2) = cds_exons(end,2)-3; 
    neg_segments = cds_exons(cds_exons(:,3)~=0,1:2);
    if genes(i).strand=='+'
      neg_segments = [neg_segments;  [cds_exons(cds_exons(:,3)==0,1)+1 cds_exons(cds_exons(:,3)==0,2)]];
      neg_segments = [neg_segments;  [cds_exons(cds_exons(:,3)==0,1)+2 cds_exons(cds_exons(:,3)==0,2)]];
      neg_segments = [neg_segments;  [cds_exons(cds_exons(:,3)==1,1)+1 cds_exons(cds_exons(:,3)==1,2)]];
      neg_segments = [neg_segments;  [cds_exons(cds_exons(:,3)==2,1)+2 cds_exons(cds_exons(:,3)==2,2)]];
    else
      neg_segments = [neg_segments;  [cds_exons(cds_exons(:,3)==0,1) cds_exons(cds_exons(:,3)==0,2)-1]];
      neg_segments = [neg_segments;  [cds_exons(cds_exons(:,3)==0,1) cds_exons(cds_exons(:,3)==0,2)-2]];
      neg_segments = [neg_segments;  [cds_exons(cds_exons(:,3)==1,1) cds_exons(cds_exons(:,3)==1,2)-1]];
      neg_segments = [neg_segments;  [cds_exons(cds_exons(:,3)==2,1) cds_exons(cds_exons(:,3)==2,2)-2]];
    end
    idx = find(neg_segments(:,2)<=neg_segments(:,1));
    neg_segments(idx,:)=[];
    idx = find(pos_segments(:,2)<=pos_segments(:,1));
    pos_segments(idx,:)=[];
  end
  pos_segments_all(num_pos+1:num_pos+size(pos_segments,1),1:size(pos_segments,2))  = pos_segments;
  neg_segments_all(num_neg+1:num_neg+size(neg_segments,1),1:size(neg_segments,2))  = neg_segments;
  
  
  num_pos = num_pos+size(pos_segments,1);
  num_neg = num_neg+size(neg_segments,1);

  last_gene_stop = gene.stop;
end

% assert(num_pos>0)
% assert(num_neg>0)

pos_segments_all(num_pos+1:end,:)=[];
neg_segments_all(num_neg+1:end,:)=[];

SEG = [pos_segments_all; neg_segments_all];
LABEL = [ones(num_pos,1); -ones(num_neg,1)];
[tmp,idx] = sort(SEG(:,1));
SEG = SEG(idx,:);
LABEL = LABEL(idx);

%assert(all(SEG(:,2)-SEG(:,1)>=0))
idx = find(SEG(:,2)-SEG(:,1)<=0);
SEG(idx,:) = [];
LABEL(idx) = [];

if ~isempty(find(SEG(2:end,1)-SEG(1:end-1,2)<0)) && ~use_frame
  idx = find(SEG(2:end,1)-SEG(1:end-1,2)<0);
  fprintf(1, 'removing %i overlapping segments\n',length(idx)*2)
  SEG(idx,:)=[];
  LABEL(idx)=[]; 
  assert(all(SEG(2:end,1)-SEG(1:end-1,2)>=0))
end
SEG = SEG(:,1:2);
