function genes = gen_gene_segmentation(P)
% genes = gen_gene_segmentation(P)
% INPUT
% P.genes
% P.segment_ids  

num_invalid_isoforms = 0;
num_isoforms = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate segmentation for all valid transcripts of a gene %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
genes = P.genes;
segment_ids = P.segment_ids;  
  
for g=1:length(genes)
  if isempty(genes(g).exons)
    continue;
  end
  if isempty(genes(g).cds_exons{1})
	genes(g) = gen_nc_gene_segmentation(genes(g), segment_ids);
	continue;
  end
  if ~genes(g).is_complete
    continue;
  end
  exons = genes(g).exons;
  for j=1:length(exons),
    if ~genes(g).transcript_valid(j)
      genes(g).segments{j}=[];
      continue;
    end
    if length(genes(g).transcript_complete)>=j&&~genes(g).transcript_complete(j)
      continue;
    end
    if genes(g).strand=='+'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % INSERT SEGMENTS FOR POSITIVE STRAND %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      isoform = [];
      % insert a row for each 5'utr exon and each intron between
      %----------------------------------------------------
      for k=1:size(genes(g).utr5_exons{j},1)
        isoform = [isoform; [genes(g).utr5_exons{j}(k,1:2) segment_ids.utr5exon]];
        if k<size(genes(g).utr5_exons{j},1)
          isoform = [isoform; [genes(g).utr5_exons{j}(k,2) genes(g).utr5_exons{j}(k+1,1) segment_ids.intron]];
        end
      end
      % if there is an intron directly before the tis site (AGATG):
      %----------------------------------------------------
      if  ~isempty(genes(g).utr5_exons{j})
        if genes(g).utr5_exons{j}(end,2)<genes(g).cds_exons{j}(1,1)
           isoform = [isoform; [genes(g).utr5_exons{j}(end,2) genes(g).cds_exons{j}(1,1) segment_ids.intron]];
        end
      end
      % insert a row for each cds_exon and each intron between
      %----------------------------------------------------
      for k=1:size(genes(g).cds_exons{j},1)
        isoform = [isoform; [genes(g).cds_exons{j}(k,1:2) segment_ids.cds_exon]];
        if k<size(genes(g).cds_exons{j},1)
          isoform = [isoform; [genes(g).cds_exons{j}(k,2) genes(g).cds_exons{j}(k+1,1) segment_ids.intron]];
        else
          % set end of last cds_exon to the beginning of the cdsStop-consensus
          isoform(end,2)=isoform(end,2)-3;
        end
      end   
      % insert an intron between last cds_exon and 3'utr if necessary
      %----------------------------------------------------   
      if  ~isempty(genes(g).utr3_exons{j}) 
        if  genes(g).cds_exons{j}(end,2)<genes(g).utr3_exons{j}(1,1)
          %isoform = [isoform; [genes(g).cds_exons{j}(end,2) genes(g).utr3_exons{j}(1,1) segment_ids.intron]];
          isoform = [isoform; [isoform(end,2) genes(g).utr3_exons{j}(1,1) segment_ids.intron]];
        end
      end
      % insert a row for each 3'utr exon and each intron between
      %----------------------------------------------------
      for k=1:size(genes(g).utr3_exons{j},1)
        if k==1 && isoform(end,3)==segment_ids.cds_exon
          % set end of first 3'utrexon to the beginning of the cdsStop-consensus
          isoform = [isoform; [genes(g).utr3_exons{j}(k,1)-3 genes(g).utr3_exons{j}(k,2) segment_ids.utr3exon]];
        % insert polya_tail if applicable
        elseif isfield(genes, 'polya') && k==size(genes(g).utr3_exons{j},1)&&length(genes(g).polya)>=j&&...
               genes(g).polya(j)>genes(g).utr3_exons{j}(k,1)&&genes(g).polya(j)<genes(g).utr3_exons{j}(k,2)
          isoform = [isoform; [genes(g).utr3_exons{j}(k,1) genes(g).polya(j) segment_ids.utr3exon]];
          isoform = [isoform; [genes(g).polya(j) genes(g).utr3_exons{j}(k,2) segment_ids.polya_tail]];
        else
          isoform = [isoform; [genes(g).utr3_exons{j}(k,1:2) segment_ids.utr3exon]];
        end
        if k<size(genes(g).utr3_exons{j},1)
          isoform = [isoform; [genes(g).utr3_exons{j}(k,2) genes(g).utr3_exons{j}(k+1,1) segment_ids.intron]];
        end
      end   
    elseif genes(g).strand=='-'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % INSERT SEGMENTS FOR NEGATIVE STRAND %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      isoform=[];
      % insert a row for each 3'utr exon and each intron between
      %----------------------------------------------------
      for k=1:size(genes(g).utr3_exons{j},1)
        % insert polya_tail if applicable
        if isfield(genes, 'polya') && k==1&&length(genes(g).polya)>=j&&genes(g).polya(j)>genes(g).utr3_exons{j}(k,1)&&...
           genes(g).polya(j)<genes(g).utr3_exons{j}(k,2)
          isoform = [isoform; [genes(g).utr3_exons{j}(k,1) genes(g).polya(j) segment_ids.polya_tail]];
          isoform = [isoform; [genes(g).polya(j) genes(g).utr3_exons{j}(k,2) segment_ids.utr3exon]];
        else
          isoform = [isoform; [genes(g).utr3_exons{j}(k,1:2) segment_ids.utr3exon]];
        end
        if k<size(genes(g).utr3_exons{j},1)
          isoform = [isoform; [genes(g).utr3_exons{j}(k,2) genes(g).utr3_exons{j}(k+1,1) segment_ids.intron]];
        else
          isoform(end,2)=isoform(end,2)+3;
        end
      end
      % if there is an intron between cdsStop and first utr3_exon:
      %----------------------------------------------------
      if  ~isempty(genes(g).utr3_exons{j})
        if genes(g).utr3_exons{j}(end,2)<genes(g).cds_exons{j}(1,1)
          isoform(end,2) = isoform(end,2)-3;%undo moving the last utr3-exon position
          isoform = [isoform; [isoform(end,2) genes(g).cds_exons{j}(1,1)+3 segment_ids.intron]];
        end
      end
      % insert a row for each cds_exon and each intron between
      %----------------------------------------------------
      for k=1:size(genes(g).cds_exons{j},1)
        if k==1 %&& ~isempty(isoform)&&isoform(end,3)==segment_ids.utr3exon
          isoform = [isoform; [genes(g).cds_exons{j}(k,1)+3 genes(g).cds_exons{j}(k,2) segment_ids.cds_exon]];
        else
          isoform = [isoform; [genes(g).cds_exons{j}(k,1:2) segment_ids.cds_exon]];
        end  
        if k<size(genes(g).cds_exons{j},1)
          isoform = [isoform; [genes(g).cds_exons{j}(k,2) genes(g).cds_exons{j}(k+1,1) segment_ids.intron]];
        end
      end  
      % insert an intron between last cds_exon and 5'utr if necessary
      %----------------------------------------------------   
      if  ~isempty(genes(g).utr5_exons{j}) 
        if genes(g).cds_exons{j}(end,2)<genes(g).utr5_exons{j}(1,1)
          isoform = [isoform; [genes(g).cds_exons{j}(end,2) genes(g).utr5_exons{j}(1,1) segment_ids.intron]];
        end
      end
      % insert a row for each 5'utr exon and each intron between
      %----------------------------------------------------
      for k=1:size(genes(g).utr5_exons{j},1)
        isoform = [isoform; [genes(g).utr5_exons{j}(k,1:2) segment_ids.utr5exon]];
        if k<size(genes(g).utr5_exons{j},1)
          isoform = [isoform; [genes(g).utr5_exons{j}(k,2) genes(g).utr5_exons{j}(k+1,1) segment_ids.intron]];
        end
      end
    end %strand case
    % assert a tight segmentation
    assert(all(isoform(2:end, 1) == isoform(1:end-1, 2)));
    if ~(all(isoform(2:end, 1) == isoform(1:end-1, 2)))
       if genes(g).utr5_exons{j}(end,2)<genes(g).tis(j)
         fprintf('last utr5_exon stop: %i \n tis: \t%i \n first cds_exon_start: %i\n\n',...
                genes(g).utr5_exons{j}(end,2),genes(g).tis(j),genes(g).cds_exons{j}(1,1))
       else 
       end
    end

    % assert that there is an ORF
    assert(sum(isoform(:,3)==segment_ids.cds_exon ...
               | isoform(:,3)==segment_ids.utr5exon ...
               | isoform(:,3)==segment_ids.utr3exon)==size(genes(g).cds_exons{j},1)+size(genes(g).utr5_exons{j},1)+size(genes(g).utr3_exons{j},1));
    cds_idx = find(isoform(:,3)==segment_ids.cds_exon);
    %assert(mod(sum(isoform(cds_idx,2)-isoform(cds_idx)),3)==0)
    if (mod(sum(isoform(cds_idx,2)-isoform(cds_idx)),3)~=0)
      %keyboard
      num_invalid_isoforms = num_invalid_isoforms+1;
    end
    num_isoforms = num_isoforms+1;
    genes(g).segments{j} = isoform;
  end
end
fprintf('%i of %i isoforms had wrong reading frame length\n', num_invalid_isoforms, num_isoforms);
return


