function  [ALL,label_2] = get_true_alternative_sites(ALL,gene,true_splice,sig,debug)

strand = gene.strand;
POS=ALL.POS;
LABEL_2=ALL.LABEL_2;

if (isequal(sig,'acc') & strand=='+')|...
      (isequal(sig,'don') & strand=='-'),
  which_end = 1;
  other_end = 2;
elseif (isequal(sig,'don')& strand=='+') | ...
      (isequal(sig,'acc')& strand=='-')
  which_end = 2;
  other_end = 1;
end
   
  
label_2 = -ones(3,length(true_splice)); 
% GET TRUE ALTERNATIVE SPLICE SITES

ii_es2 = []; ii_o2 = [];  ii_ir2 = [];  ii_3p2 = [];  ii_5p2 = []; 

if gene.alt_exon %% EXON SKIP
  alt_exons = gene.alt_exon_details.exons;
  true_alt_es = gene.splicegraph{1}(which_end,alt_exons(2,:)) ;
  [temp,ii1,ii_es2] = intersect(true_alt_es,true_splice);
  assert(length(temp)==length(unique(true_alt_es)))
  true_alt_es = true_alt_es(ii1);
  % splice_before_es = gene.splicegraph{1}(don_end,alt_exons(1,ii1)) ;
  label_2(1,ii_es2) = 1;
  if debug
    fprintf(1,'exon skip')
    viewsplicegraph(gene,1)
    hold on
    for jj=1:length(true_alt_es)
      plot([true_alt_es(jj) true_alt_es(jj)],[0 10],'b');
    end
    if keyboard_allowed(),
      keyboard
    end ;
  end
end %if gene.alt_exon

if gene.alt_intron
  alt_introns = gene.alt_intron_details.exons;
  true_alt_ir = gene.splicegraph{1}(which_end,alt_introns(other_end,:));
  [temp,ii1,ii_ir2] = intersect(true_alt_ir,true_splice);
  assert(length(temp)==length(unique(true_alt_ir)))
  true_alt_ir = true_alt_ir(ii1);
  % splice_before_ir = gene.splicegraph{1}(acc_end,alt_introns(1,ii1));
  % intron_ret(1,:) = gene.splicegraph{1}(which_end,alt_introns(other_end,ii1));
  % intron_ret(2,:) = gene.splicegraph{1}(other_end,alt_introns(which_end,ii1));
  label_2(2,ii_ir2) = 1;
  if debug
    fprintf(1,'intron retention')
    viewsplicegraph(gene,1)
    hold on
    for jj=1:length(true_alt_ir)
      plot([true_alt_ir(jj) true_alt_ir(jj)],[0 10],'b');
    end
    if keyboard_allowed(),
      keyboard
    end 
  end
end

if gene.alt_5prime&isequal(sig, 'don') 
  %if ~isempty(gene.alt_5prime_details), keyboard; end ;
  if gene.strand=='+',
    alt_5p = gene.alt_5prime_details.exons(1,:);
  else
    alt_5p = gene.alt_5prime_details.exons(2,:);
  end ;
  alt_5p = unique([alt_5p gene.alt_5prime_details.exons(3,:)]);
  true_alt_5p = gene.splicegraph{1}(which_end,alt_5p);
  [temp,ii1,ii_5p2] = intersect(true_alt_5p,true_splice);
  assert(length(temp)==length(unique(true_alt_5p)))
  true_alt_5p = true_alt_5p(ii1);
  % splice_before_5p = gene.splicegraph{1}(acc_end,);
  % don_after = gene.splicegraph{1}(other_end,alt_3p);
  label_2(3,ii_5p2) = 1;
  if debug
    fprintf(1,'alt 5prime')
    viewsplicegraph(gene,1)
    hold on
    for jj=1:length(true_alt_5p)
      plot([true_alt_5p(jj) true_alt_5p(jj)],[0 10],'b');
    end
    if keyboard_allowed(),
      keyboard
    end 
  end
end
if gene.alt_3prime&isequal(sig,'acc')
  %if ~isempty(gene.alt_3prime_details), keyboard; end ;
  if gene.strand=='+',
    alt_3p = gene.alt_3prime_details.exons(2,:);
  else
    alt_3p = gene.alt_3prime_details.exons(1,:);
  end ;
  alt_3p = unique([alt_3p gene.alt_3prime_details.exons(3,:)]);
  true_alt_3p = gene.splicegraph{1}(which_end,alt_3p);
  [temp,ii1,ii_3p2] = intersect(true_alt_3p,true_splice);
  assert(length(temp)==length(unique(true_alt_3p)))
  true_alt_3p = true_alt_3p(ii1);
  % splice_before_3p = gene.splicegraph{1}(don_end,[exon_alt_3prime(find(idx_alt_3prime==g_idx)).fiveprimesite]);
  % don_before = gene.splicegraph{1}(other_end,[exon_alt_3prime(find(idx_alt_3prime==g_idx)).fiveprimesite]);
  % don_after = gene.splicegraph{1}(other_end,alt_3p);
  label_2(3,ii_3p2) = 1;
  if debug
    fprintf(1,'alt 3prime')
    viewsplicegraph(gene,1)
    hold on
    for jj=1:length(true_alt_3p)
      plot([true_alt_3p(jj) true_alt_3p(jj)],[0 10],'b');
    end
    if keyboard_allowed(),
      keyboard 
    end 
  end
end
[ii,ii1,ii2] = intersect(POS,true_splice);
LABEL_2(:,ii1) = label_2(:,ii2);




