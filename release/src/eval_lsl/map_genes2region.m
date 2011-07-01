function [nuc,exons,transcripts,utr5, utr3,signals] = map_genes2region(genes,region,mask_regions,filter_utr)
%  [nuc,exons,transcripts] = map_genes2region(genes,region,mask_regions)
  
nuc.all  = zeros(1,region.stop-region.start+1) ;
nuc.cds = zeros(1,region.stop-region.start+1) ;

exons.all = zeros(length(genes)*10,2);
exons.cds =  zeros(length(genes)*10,2);
exons.all_gene_ids =  zeros(length(genes)*10,1); 
exons.cds_gene_ids = zeros(length(genes)*10,1);
num.exons = 0;
num.cds_exons = 0;

transcripts.all={};
transcripts.cds={};
transcripts.status_all=[];
transcripts.status_cds=[];

utr5 = zeros(length(genes)*10,2);
utr3 = zeros(length(genes)*10,2);
num.utr5p =0;
num.utr3p =0;

sig_names = {'tis','cdsStop','tss','cleave','acc','don'};
for i=1:length(sig_names)
  signals.(sig_names{i}) = zeros(length(genes)*10,1);
  num.(sig_names{i})=0;
end

gene_id  = [];
cds_gene_id  = [];
num_all = 0;
num_cds = 0;

strand = unique([genes.strand]);
assert(length(strand)==1 || isempty(genes))


for g_idx=1:length(genes)
  gene = genes(g_idx);
  if mod(g_idx,100)==0
    fprintf('region %i / %i \r',g_idx,length(genes))
  end
  %% signals
  if isfield(gene, 'tis') && isfield(gene, 'tss') && isfield(gene, 'cleave') && isfield(gene, 'cdsStop') && iscell(gene.tis)
    signals.tis(num.tis+[1:length([gene.tis{:}])]) = [gene.tis{:}] ;
    signals.cdsStop(num.cdsStop+[1:length([gene.cdsStop{:}])]) = [gene.cdsStop{:}] ;
    num.tis = num.tis + length([gene.tis{:}]);
    num.cdsStop = num.cdsStop + length([gene.cdsStop{:}]);
    tss = [];
    for j=1:length(gene.tis)
      if ~isempty(gene.tis{j}) && length(gene.tss)>=j && ~isempty(gene.tss{j}) && abs(gene.tis{j}-gene.tss{j})>=filter_utr
        tss = [tss gene.tss{j}];
        utr5(num.utr5p+1,:) = [gene.tss{j} gene.tis{j}];
        num.utr5p = num.utr5p+1;
      end
    end
    cleave = [];
    for j=1:length(gene.cdsStop)
      if ~isempty(gene.cdsStop{j}) && length(gene.cleave)>=j && ~isempty(gene.cleave{j}) && abs(gene.cdsStop{j}-gene.cleave{j})>=filter_utr;    
        cleave = [cleave gene.cleave{j}];
        utr3(num.utr3p+1,:) = [gene.cdsStop{j} gene.cleave{j}];
        num.utr3p = num.utr3p+1;
      end
    end
    signals.tss(num.tss+[1:length(tss)]) = tss;
    signals.cleave(num.cleave+[1:length(cleave)]) = cleave;
    num.tss = num.tss + length(tss);
    num.cleave = num.cleave + length(cleave);
  end
 
  for j=1:length(gene.exons)
    %exons
    %exons.all(num.exons+[1:size(gene.exons{j},1)],1:2) = gene.exons{j}(:,1:2);
    %exons.all(num.exons+[1:size(gene.exons{j},1)-2],1:2) = gene.exons{j}(2:end-1,1:2);
    tmp=gene.exons{j}(2:end-1,1:2) ;
    if ~isempty(tmp),
      if ~(gene.exons{j}(1,1)<gene.exons{j}(2,1)) 
        warning('exons are misordered')
      end
      if ~(gene.exons{j}(1,1)<gene.exons{j}(1,2))
        warning('exons are misordered')
      end
      tmp=[gene.exons{j}(1,2)-10 gene.exons{j}(1,2); tmp; gene.exons{j}(end,1) gene.exons{j}(end,1)+10] ;
    end ;
    exons.all(num.exons+[1:size(tmp,1)],1:2) = tmp;
    exons.all_gene_ids(num.exons+[1:size(gene.exons{j},1)])  = gene.id ;
    exons.cds(num.cds_exons+[1:size(gene.cds_exons{j},1)],1:2) = gene.cds_exons{j}(:,1:2);
    exons.cds_gene_ids(num.cds_exons+[1:size(gene.cds_exons{j},1)])  = gene.id;
    
    num.exons = num.exons +size(gene.exons{j},1);
    num.cds_exons = num.cds_exons +size(gene.cds_exons{j},1);
    len = length(signals.acc);
    %% splice sites
    if strand =='+'
      signals.acc(num.acc+(1:size(gene.cds_exons{j},1)-1)) = gene.cds_exons{j}(2:end,1);
      signals.don(num.don+(1:size(gene.cds_exons{j},1)-1)) = gene.cds_exons{j}(1:end-1,2);
    else
      signals.acc(num.acc+(1:size(gene.cds_exons{j},1)-1)) = gene.cds_exons{j}(1:end-1,2);
      signals.don(num.don+(1:size(gene.cds_exons{j},1)-1)) = gene.cds_exons{j}(2:end,1);
    end
    %if length(signals.acc)>len
    %  warning('array was initialized to small: loosing time copying everything\n');
    %  keyboard
    %end
	if length(signals.acc)>len
		signals.acc = [signals.acc; ones(len, 1)];
		signals.don = [signals.don; ones(len, 1)];
	end

    num.acc = num.acc + max(size(gene.cds_exons{j},1)-1, 0) ;
    num.don = num.don + max(size(gene.cds_exons{j},1)-1, 0) ;

    % transcripts (all exons)
    num_all = num_all + 1;
    transcripts.all{num_all}=gene.exons{j}(:,1:2);
	
    % transcripts (cds exons)
    num_cds = num_cds + 1;
    transcripts.cds{num_cds}=gene.cds_exons{j}(:,1:2);

    if isfield(gene, 'transcript_status')
      if ~isempty(gene.transcript_status)&length(gene.transcript_status)>=j
        transcripts.status_all(num_all) = gene.transcript_status(j);
      elseif~isempty(gene.transcript_status)
        transcripts.status_all(num_all) = gene.transcript_status(1);   
      end

      if ~isempty(gene.transcript_status)&length(gene.transcript_status)>=j
        transcripts.status_cds(num_cds) = gene.transcript_status(j);
      elseif~isempty(gene.transcript_status)
        transcripts.status_cds(num_cds) = gene.transcript_status(1);   
      end
    end
  end%loop over transcripts
end

all_unique_idx = unique_transcripts(transcripts.all);
cds_unique_idx = unique_transcripts(transcripts.cds);

transcripts.all = transcripts.all(find(all_unique_idx==1))';
transcripts.cds = transcripts.cds(find(cds_unique_idx==1))';


for i=1:length(sig_names)
  signals.(sig_names{i})(num.(sig_names{i})+1:end) = [];  
  signals.(sig_names{i}) = unique(signals.(sig_names{i}));
end
exons.all(num.exons+1:end,:) = [];  
exons.cds(num.cds_exons+1:end,:) = [];  
exons.all_gene_ids(num.exons+1:end) = [];  
exons.cds_gene_ids(num.cds_exons+1:end) = [];  

[exons.all,idx] = unique(exons.all,'rows');
exons.all_gene_ids= exons.all_gene_ids(idx);
[exons.cds,idx] = unique(exons.cds,'rows');
exons.cds_gene_ids = exons.cds_gene_ids(idx);

utr3(num.utr3p+1:end,:) = [];  
utr5(num.utr5p+1:end,:) = [];  

utr3 = unique(utr3,'rows');
utr5 = unique(utr5,'rows');

if ~isempty(mask_regions)
  rm_exon=[];
  for k=1:size(exons.all,1)
    if any(mask_regions(1,:)<=exons.all(k,1)&mask_regions(2,:)>=exons.all(k,2))%remove exons that are included in mask_region 
      rm_exon=[rm_exon k];
    elseif any(mask_regions(1,:)>=exons.all(k,1)&mask_regions(1,:)<=exons.all(k,2)) %shorten exons that overlap mask_region 
      idx = find(mask_regions(1,:)>=exons.all(k,1)&mask_regions(1,:)<=exons.all(k,2));
      exons.all(k,2) = min(mask_regions(1,idx));
    elseif any(mask_regions(2,:)>=exons.all(k,1)&mask_regions(2,:)<=exons.all(k,2))
      idx = find(mask_regions(2,:)>=exons.all(k,1)&mask_regions(2,:)<=exons.all(k,2));
      exons.all(k,1) = max(mask_regions(2,idx));
    end
  end
  exons.all(rm_exon,:)=[];
  rm_exon=[];
  for k=1:size(exons.cds,1)
    if any(mask_regions(1,:)<=exons.cds(k,1)&mask_regions(2,:)>=exons.cds(k,2))%remove exons that are included in mask_region 
      rm_exon=[rm_exon k];
    elseif any(mask_regions(1,:)>=exons.cds(k,1)&mask_regions(1,:)<=exons.cds(k,2)) %shorten exons that overlap mask_region 
      idx = find(mask_regions(1,:)>=exons.cds(k,1)&mask_regions(1,:)<=exons.cds(k,2));
      exons.cds(k,2) = min(mask_regions(1,idx));
    elseif any(mask_regions(2,:)>=exons.cds(k,1)&mask_regions(2,:)<=exons.cds(k,2))
      idx = find(mask_regions(2,:)>=exons.cds(k,1)&mask_regions(2,:)<=exons.cds(k,2));
      exons.cds(k,1) = max(mask_regions(2,idx));
    end
  end
  exons.cds(rm_exon,:)=[];
  warning('only implemented for exon level')
end


for k=1:size(exons.all,1)
  start = max(1,exons.all(k,1)-region.start+1);
  stop = min(length(nuc.all),exons.all(k,2)-region.start+1);
  nuc.all([start:stop])=1;
end
for k=1:size(exons.cds,1)
  start = max(1,exons.cds(k,1)-region.start+1);                                                                        
  stop = min(length(nuc.cds),exons.cds(k,2)-region.start+1);  
  nuc.cds([start:stop])=1;
end


nuc.all(nuc.all~=0)=1;
nuc.cds(nuc.cds~=0)=1;

nuc.all = find(nuc.all==1)';
nuc.cds = find(nuc.cds==1)';
