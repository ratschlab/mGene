function genes = cover_splicegraph(genes)
% function genes = cover_splicegraph(genes)
% 
% From the results of annotate_genes, "cover" the splicegraph with
% known alternative splicing and transcription events.
% 
% Written by: Cheng Soon Ong, 31 May 2007

for ix = 1:length(genes)
  if mod(ix,1000)==0
    fprintf('.');
  end

  % 0 means not used, 1 means used, > 1 means alternative.
  isaltv = ones(1,size(genes(ix).splicegraph{1},2));
  isalte = genes(ix).splicegraph{2};
  vertices = genes(ix).splicegraph{1};
  edges = genes(ix).splicegraph{2};

  if genes(ix).alt_exon > 0
    for ixe = 1:size(genes(ix).alt_exon_details.exons,2)
      exons = [genes(ix).alt_exon_details.exons(:,ixe)];
      assert(edges(exons(1),exons(2))==1 && edges(exons(2),exons(3))==1 && edges(exons(1),exons(3))==1);
      isaltv(exons) = isaltv(exons) + 1;
      isalte(exons(1),exons(2)) = isalte(exons(1),exons(2)) + 1;
      isalte(exons(2),exons(1)) = isalte(exons(2),exons(1)) + 1;
      isalte(exons(2),exons(3)) = isalte(exons(2),exons(3)) + 1;
      isalte(exons(3),exons(2)) = isalte(exons(3),exons(2)) + 1;
      isalte(exons(1),exons(3)) = isalte(exons(1),exons(3)) + 1;
      isalte(exons(3),exons(1)) = isalte(exons(3),exons(1)) + 1;
    end
  end

  if genes(ix).alt_intron > 0
    for ixe = 1:size(genes(ix).alt_intron_details.exons,2)
      exons = [genes(ix).alt_intron_details.exons(:,ixe)];
      assert(edges(exons(1),exons(2))==1)
      isaltv(exons) = isaltv(exons) + 1;
      isalte(exons(1),exons(2)) = isalte(exons(1),exons(2)) + 1;
      isalte(exons(2),exons(1)) = isalte(exons(2),exons(1)) + 1;
      tedges= triu(edges) ;
      if all(tedges(:,exons(1))==tedges(:,exons(3))),
        isalte(tedges(:,exons(1))~=0,exons(1)) = isalte(tedges(:,exons(1))~=0,exons(1)) + 1 ;
        isalte(tedges(:,exons(1))~=0,exons(3)) = isalte(tedges(:,exons(1))~=0,exons(3)) + 1 ;
        isalte(exons(1),tedges(:,exons(1))~=0) = isalte(exons(1),tedges(:,exons(1))~=0) + 1 ;
        isalte(exons(3),tedges(:,exons(1))~=0) = isalte(exons(3),tedges(:,exons(1))~=0) + 1 ;
      end ;
      tedges= tril(edges) ;
      if all(tedges(:,exons(2))==tedges(:,exons(3))),
        isalte(tedges(:,exons(2))~=0,exons(1)) = isalte(tedges(:,exons(2))~=0,exons(2)) + 1 ;
        isalte(tedges(:,exons(2))~=0,exons(3)) = isalte(tedges(:,exons(2))~=0,exons(3)) + 1 ;
        isalte(exons(2),tedges(:,exons(2))~=0) = isalte(exons(2),tedges(:,exons(2))~=0) + 1 ;
        isalte(exons(3),tedges(:,exons(2))~=0) = isalte(exons(3),tedges(:,exons(2))~=0) + 1 ;
      end ;
      
    end
  end

  if genes(ix).alt_5prime > 0
    for ixe = 1:size(genes(ix).alt_5prime_details.exons,2)
      exons = [genes(ix).alt_5prime_details.exons(:,ixe)];
      assert(edges(exons(1),exons(2))==1);
      isaltv(exons) = isaltv(exons) + 1;
      isalte(exons(1),exons(2)) = isalte(exons(1),exons(2)) + 1;
      isalte(exons(2),exons(1)) = isalte(exons(2),exons(1)) + 1;
      if genes(ix).strands(1) == '-'
        assert(edges(exons(1),exons(3))==1);
        isalte(exons(1),exons(3)) = isalte(exons(1),exons(3)) + 1;
        isalte(exons(3),exons(1)) = isalte(exons(3),exons(1)) + 1;
      else
        assert(edges(exons(2),exons(3))==1);
        isalte(exons(2),exons(3)) = isalte(exons(2),exons(3)) + 1;
        isalte(exons(3),exons(2)) = isalte(exons(3),exons(2)) + 1;
      end
    end
  end

  if genes(ix).alt_3prime > 0
    for ixe = 1:size(genes(ix).alt_3prime_details.exons,2)
      exons = [genes(ix).alt_3prime_details.exons(:,ixe)];
      assert(edges(exons(1),exons(2))==1);
      isaltv(exons) = isaltv(exons) + 1;
      isalte(exons(1),exons(2)) = isalte(exons(1),exons(2)) + 1;
      isalte(exons(2),exons(1)) = isalte(exons(2),exons(1)) + 1;
      if genes(ix).strands(1) == '+'
	assert(edges(exons(1),exons(3))==1);
	isalte(exons(1),exons(3)) = isalte(exons(1),exons(3)) + 1;
	isalte(exons(3),exons(1)) = isalte(exons(3),exons(1)) + 1;
      else
	assert(edges(exons(2),exons(3))==1);
	isalte(exons(2),exons(3)) = isalte(exons(2),exons(3)) + 1;
	isalte(exons(3),exons(2)) = isalte(exons(3),exons(2)) + 1;
      end
    end
  end

  % if region is alternative, but isaltv and isalte doesn't think so,
  % call region complex.
  num_alt_complex = 0;
  gene = genes(ix) ;
  if isempty(genes(ix).alt_tstart_details),
    genes(ix).alt_tstart_details.exons=[] ;
  end ;
  if isempty(genes(ix).alt_tend_details),
    genes(ix).alt_tend_details.exons=[] ;
  end ;
  idx=setdiff(1:size(gene.splicegraph{1},2), union(genes(ix).alt_tstart_details.exons,genes(ix).alt_tend_details.exons)) ;
  gene.splicegraph{1}=gene.splicegraph{1}(:,idx) ;
  gene.splicegraph{2}=gene.splicegraph{2}(idx,idx) ;
  vertices = gene.splicegraph{1};
  edges = gene.splicegraph{2};
  region = detectaltregions(gene);
  isaltv_ = isaltv(idx) ;
  isalte_ = isalte(idx,idx) ;
  for ixr = 1:size(region,2)
    assert(region(3,ixr)>1);
    exons = find([vertices(1,:)<region(2,ixr)] & [vertices(2,:)>region(1,ixr)]);
    if any([isaltv_(exons) == 1]) || any(any([isalte_(exons,exons) == 1]))
      %keyboard;
      num_alt_complex = num_alt_complex + 1;
      genes(ix).exons_is_alt(exons) = genes(ix).exons_is_alt(exons) + 1;
    end
  end
  genes(ix).num_alt_complex = num_alt_complex;
end

%
if 0,
  
  if genes(ix).alt_tstart > 0
    for exons = genes(ix).alt_tstart_details.exons,
      isaltv(exons) = isaltv(exons) + 1;
      isalte(exons,find(edges(exons,:))) = isalte(exons,find(edges(exons,:))) + 1 ;
      isalte(find(edges(:,exons)),exons) = isalte(find(edges(:,exons)),exons) + 1 ;
      for next_exons = find(edges(exons,:))
        istart = vertices(2,exons) ;
        iend = vertices(1,next_exons) ;
        for i=1:size(edges,1),
          for j=i+1:size(edges,2),
            if vertices(2,i)==istart && vertices(1,j)==iend,
              isalte(i,j) = isalte(i,j)+1 ;
              isalte(j,i) = isalte(j,i)+1 ;
            end ;
          end ;
        end ;
      end ;
    end ;
  end

  if genes(ix).alt_tend > 0
    for exons = genes(ix).alt_tend_details.exons,
      isaltv(exons) = isaltv(exons) + 1;
      isalte(exons,find(edges(exons,:))) = isalte(exons,find(edges(exons,:))) + 1 ;
      isalte(find(edges(:,exons)),exons) = isalte(find(edges(:,exons)),exons) + 1 ;
    end ;
  end

end ;