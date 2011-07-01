function [ALL,is_invalid] = basic_cand_stuff(gene,ALL,sig,debug,signal_name)
% [ALL,is_invalid] = basic_cand_stuff(gene,ALL,sig,debug,signal_name)

if nargin<5
  signal_name =[];
end

is_invalid=0; 

genic_idx = find(ALL.POS>gene.start&ALL.POS<gene.stop);
ALL.GENE_ID(genic_idx) = gene.id;

if isfield(gene,'transacc') && ~isempty(gene.transacc)
  ALL.TRANS(genic_idx) = 1;
end
if ~isempty(gene.tss)
  ALL.TSS(genic_idx) = 1;
end

if ~isempty(signal_name)&&~isempty(gene.(signal_name))
  ALL.HAS_SIG(genic_idx) = 1;
end
  
if ~gene.is_valid
  is_invalid=1; 
  ALL.GENE_ID(genic_idx) =  - gene.id;
  return
end

if isequal(sig,'other')
  if gene.is_alt 
    ALL.ALTGENIC(genic_idx) = 1;
  end
else
  if any(gene.exon_in_alt_region)
    alt_reg = get_connected_regions(gene.exon_in_alt_region);  
    
    for i=1:length(unique(alt_reg))-1
      temp=find(alt_reg==i);
      alt_start = gene.splicegraph{1}(1,temp(1));
      alt_stop = gene.splicegraph{1}(2,temp(end));
      alt_genic_idx = find(ALL.POS>alt_start&ALL.POS<alt_stop);
      ALL.ALTGENIC(alt_genic_idx) = 1;
      if debug
        viewsplicegraph(gene,1)
        hold on 
        plot([alt_start alt_start],[1 10],'b')
        plot([alt_stop alt_stop],[1 10],'b')  
      end
    end
  end
end

if gene.is_alt 
  ALL.ALTGENIC(genic_idx) = 1;
end