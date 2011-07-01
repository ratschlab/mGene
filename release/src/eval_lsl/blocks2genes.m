function [genes,operons] = blocks2genes(blocks,model,name,genome_info,clade)
% [genes,operons] = blocks2genes(blocks,model,name,genome_info,clade)
%name='prediction' , name='truth'

cnt_genes = 0;
PAR=init_PAR;
if nargin>4 & isequal(clade,'nematode')
  PAR.organism.clade='nematode';
else
  clade=[];
end
% genes = init_genes(PAR);
operons = init_regions; 
tail = [];utr3=[];
cnt_operons = 0;
op=[];

intergenic_seg_id = model.segments.intergenic;
if isfield(model.segments,'intergenictrans') 
  intergenic_seg_id = [intergenic_seg_id model.segments.intergenictrans] ;
end
if isfield(model.segments, 'intercistronic')
  intergenic_seg_id = [intergenic_seg_id model.segments.intercistronic] ;
end

for i = 1:length(blocks)
  block = blocks(i);
  if isempty(block.(name))
    fprintf('block %i',i)
    warning('block empty')
    continue
  end
  
  if ~isfield(block.(name),'genes')|isequal(clade,'nematode')
    if ~isfield(block.(name),'segments')
      [block.(name).segments block.(name).genes]  = path2segmentation(block.(name).path, block.(name).pos_idx, block, model);
    end
    segments = block.(name).segments;

    idx =  find(ismember(segments(:,3),intergenic_seg_id));
    n_genes = length(idx)-1;
    block.(name).genes{1} = [];
    internal_gene = zeros(1,n_genes+1);
    cnt = 0;
    for j=1:n_genes
      if idx(j)+2==idx(j+1)&isinf(segments(idx(j)+1,3));
        continue
      elseif isinf(segments(idx(j)+1,3));
        keyboard
      end
      cnt=cnt+1;
      block.(name).genes{cnt} = segments(idx(j)+1:idx(j+1)-1,:);
      internal_gene(cnt)=segments(idx(j),3)==model.segments.intercistronic;
    end
    n_genes = length(idx)-1-sum(isinf(segments(:,3)));
    % assert(n_genes==length(block.(name).genes) | (n_genes==0 && length(block.(name).genes)==1 && isempty(block.(name).genes{1})))
    op = [];
    operon_starts = find(internal_gene(2:end)-internal_gene(1:end-1)==1);
    operon_stops = find(internal_gene(1:end-1)-internal_gene(2:end)==1);
    n_operons =  length(operon_starts);
    cnt_operons = cnt_operons+n_operons;
    for j=1:n_operons
      op(j).genes = [operon_starts(j):operon_stops(j)];
    end
  end

  
  if isempty(block.(name).genes{1})
    continue
  end
  all_transcripts = {};
  for j=1:length(block.(name))
    all_transcripts=[all_transcripts{:} block.(name)(j).genes];
  end
  mRNA = init_genes(length(all_transcripts),PAR);
  
  % fprintf('number of transcripts in block %i (reg %i%s):\t %i\n',i,block.reg,block.strand,length(mRNA));
  
  for j = 1:length(all_transcripts)
    exons =  all_transcripts{j};
    utr5_exon = exons(find(exons(:,3)==model.segments.utr5exon),1:3) ;
    idx=find(utr5_exon(:,2)==utr5_exon(:,1));
    if ~isempty(idx)
      warning('AGATG case!!')
      utr5_exon(idx,:) = [];
    end
   
    if isfield(model.segments,'transexon')
      trans_exon = exons(find(exons(:,3)==model.segments.transexon),1:3) ;
      % trans_exon(:,3) = model.segments.utr5exon;
    end
    cds_exon = exons(find(exons(:,3)==model.segments.cds_exon),1:3) ;
    %if isempty(cds_exon)
    %  cds_exon = [exons(find(isnan(exons(:,3))),1:2) model.segments.cds_exon] ;
    %end ;
    utr3_exon = exons(find(exons(:,3)==model.segments.utr3exon),1:3) ;
    if isempty(utr3_exon)
      if ~isempty(cds_exon),
      	warning('empty utr3_exon');
        utr3_exon = [cds_exon(end,2) cds_exon(end,2) model.segments.utr3exon];
      end ;
    end
	if ~isempty(utr3_exon)
    	utr3 = [utr3,max(max(utr3_exon(:,1:2)))-min(min(utr3_exon(:,1:2)))];
	end
    polya_tail = exons(find(exons(:,3)==model.segments.polya_tail),1:3) ;
    if ~isempty(polya_tail)
      assert(size(polya_tail,1)==1)
      tail = [tail,polya_tail(1,2)-polya_tail(1,1)];
      idx=find(utr3_exon(:,2)==polya_tail(1,1));
      utr3_exon(idx,2)=polya_tail(1,2);
    end
    [tmp,idx] = max(cds_exon(:,2));
    idx2 = find(utr3_exon(:,1)==cds_exon(idx,2));
    utr3_exon(idx2,1) = utr3_exon(idx2,1)+3;
    if utr3_exon(idx2,1)>=utr3_exon(idx2,2)
      utr3_exon(idx2,2)=utr3_exon(idx2,1)+1;
    end
    cds_exon(idx,2) = cds_exon(idx,2)+3;
        
    clear exons
    % keyboard
    if isfield(model.segments,'transexon')
      exons = [trans_exon; utr5_exon ;cds_exon ;utr3_exon] ;
    else
      exons = [utr5_exon ;cds_exon ;utr3_exon] ;
    end
	if isempty(exons) && isfield(model.use, 'non_coding') && model.use.non_coding
		exons =  all_transcripts{j};
		exons = exons(find(exons(:,3)==model.segments.nc_exon),1:3) ;
	end
    %%     transform to genome coordinates
    if block.strand =='+'
      exons(:,1:2)= exons(:,1:2)+block.offset ;
      polya_tail(:,1:2) = polya_tail(:,1:2) + block.offset ;
      if ~isempty(polya_tail)
        mRNA(j).polya =  polya_tail(1,1);
      end
    else
      exons(:,1:2) = -exons(:,1:2)+block.offset ;
      exons = exons(end:-1:1,:) ;
      exons(:,1:2) =exons(:,2:-1:1);
      polya_tail(:,1:2) = -polya_tail(:,1:2) + block.offset ;
      polya_tail = polya_tail(end:-1:1,:) ;
      polya_tail(:,1:2) =polya_tail(:,2:-1:1);
      if ~isempty(polya_tail)
        mRNA(j).polya =  polya_tail(1,2);
      end
    end
    if isfield(model.segments,'transexon')
      idx = find(exons(:,3)==model.segments.transexon);
      assert(isempty(idx)|length(idx)==1);
      mRNA(j).transacc = exons(idx,1);
      %       exons(idx,:)= [];
    end
    
    mRNA(j).exons{1} = exons;
    mRNA(j).start = min(min(exons(:,1:2))); 
    mRNA(j).stop = max(max(exons(:,1:2)));

    mRNA(j).strand = block.strand;
    mRNA(j).chr = block.chr;
    if isfield(block,'chr_num')&~isempty(block.chr_num)
      mRNA(j).chr_num = block.chr_num;
    end
    tname = sprintf('%s%s_gene=%i',block.chr,block.strand,cnt_genes+j);
    mRNA(j).name = tname;
    mRNA(j).transcripts{1} = tname;
    
    if isequal(clade,'nematode')& ~isempty(op)&& ismember(j,[op.genes])
      for k=1:length(op)
        if ismember(j,op(k).genes)
          mRNA(j).in_operon = find(op(k).genes==j);
        end
      end
    elseif isequal(clade,'nematode') & ~isempty(op)&&  ~ismember(j,[op.genes])
      mRNA(j).in_operon = 0;
    end
    
  end
  genes(cnt_genes+[1:length(mRNA)])=mRNA;
  if ~isempty(op)
    for j=1:length(op)
      operons(end+1).chr = block.chr;
      if isfield(block,'chr_num')&~isempty(block.chr_num)
        operons(end).chr_num = block.chr_num;
      end
      operons(end).strand = block.strand;
      operons(end).start =  min([mRNA(op(j).genes).start]);
      operons(end).stop =  max([mRNA(op(j).genes).stop]);
      operons(end).id = length(operons)-1;
      operons(end).num_genes = length(op(j).genes);
    end
  end
  
  cnt_genes = cnt_genes+length(mRNA);
  fprintf('number of genes %i\n',cnt_genes)

end
operons = operons(2:end);

if isfield(model.use, 'non_coding') && model.use.non_coding
	coding = zeros(1, length(genes));
	for j=1:length(genes)
		if ismember(model.segments.cds_exon, genes(j).exons{1}(:, 3))
			coding(j) = 1;
		end
	end
	coding_idx = find(coding);
	genes(coding_idx) = process_gff3_genes(genes(coding_idx), genome_info, 'CDS', model) ; 
else
	genes = process_gff3_genes(genes, genome_info, 'CDS', model) ;
end
which_checks.exons_sorted = 1;
which_checks.intron_length = 1;
which_checks.splicesites = 1;
which_checks.orf = 1;
which_checks.gene_length = 1;
which_checks.graph = 0;
which_checks.transacc = 0;

%genes = check_genes(genes, genome_info,which_checks);

% build splicegraphs
%genes = build_splice_graph_caller(genes);
%genes = infer_splice_graph_caller(genes) ;

%fprintf('Sorting exons...\n');
%for i=1:length(genes)
%  if isempty(genes(i).splicegraph{1}), 
%    continue ; 
%  end ;
%  [dummy,exon_order] = sort(genes(i).splicegraph{1}(1,:),2,'ascend');
%  genes(i).splicegraph{1} = genes(i).splicegraph{1}(:,exon_order);
%  genes(i).splicegraph{2} = genes(i).splicegraph{2}(exon_order,exon_order);
%end
%genes = alt_const(genes);
%
%length(tail)
% hist(tail,100)
