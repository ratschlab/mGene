function genes = correct_tis_stop(genes, genome_info, min_orf_len, min_orf_sep, verb) ;
% genes = correct_tis_stop(genes, genome_info, min_orf_len, min_orf_sep, verb) ;

if nargin==1,
  PAR=genes ;

  genes=PAR.genes ;
  genome_info = PAR.genome_info ;
  min_orf_len = PAR.min_orf_len ;
  min_orf_sep = PAR.min_orf_sep ;
  verb = PAR.verb ;

else
  if nargin<3,
    min_orf_len = 300 ;
  end ;
  if nargin<4,
    min_orf_sep = 0.7 ;
  end ;
  if nargin<5,
    verb=2 ;
  end ;
end ;

num_tis_shift_max = 0 ;
num_stop_shift_max = 0 ;

length(genes)
for i=1:length(genes), 
  if verb==1 && mod(i,100)==0,
    fprintf('%i/%i ... ', i, length(genes)) ;
  end ;

  if isempty(genes(i).exons)||all(cellfun('isempty',genes(i).exons))
    continue
  end
  if isempty(genes(i).transcripts)||all(cellfun('isempty',genes(i).transcripts))
    continue
  end
  
  %%% CHECK ORF
  for j=1:length(genes(i).exons), 

    exons = genes(i).exons{j};
    if isempty(exons) || exons(1,1)==0 ;
      continue ;
    end ;
    if ~isempty(genes(i).cds_exons{j}),
      continue ;
    end ;

    if genes(i).strand=='+'
      exons(:,2) = exons(:,2)-1;
    else
      exons(:,1) = exons(:,1)+1;
    end
    [genestr,tmp,ok] = load_genomic(genes(i).chr, genes(i).strand, exons(:,1), exons(:,2), genome_info, 1);
    %if any(any(tmp(2:end-1,:)==0)),
    %  fprintf('gene %i has some problematic splice sizes\n', i) ;
    %  %keyboard
    %end ;
    [s, e, atg]=find_max_orfs(genestr, 'atg') ;
    cand = find(atg~=0 & e<length(genestr)-2) ; % exclude cases with no atg or stop
    cand2 = find(atg~=0) ; % exclude cases with no atg or stop
    [cds_len, idx] = max(e(cand)-atg(cand)) ;
    if length(cds_len)~=1 
      if verb>=2,
        fprintf('=> Cannot suitable identify ORF\n') ;
      end ;
      continue ;
    end ;
    cds_lens = sort(e(cand2)-atg(cand2)) ;
    if length(cds_lens)>1 && cds_len*min_orf_sep<max(cds_lens(end-1)),
      if verb>=2,
        fprintf('=> Cannot unambiguously identify ORF (%i vs %i)\n', cds_len, cds_lens(end-1)) ;
      end ;
      continue ;
    end ;
    
    if ~isempty(cand) && ~isempty(idx) && cds_len>min_orf_len,


      tis_pos     = map_rna_pos(genes(i).exons{j}, genes(i).strand, atg(cand(idx)), 0) ;
      cdsStop_pos = map_rna_pos(genes(i).exons{j}, genes(i).strand, e(cand(idx))+3+1, 1) ;
      
      %keyboard
      
      [cds_exons, utr5_exons, utr3_exons, ok] = split_exons_cds(genes(i).exons{j}, tis_pos, cdsStop_pos, genes(i).strand) ;
      if ok,
        assert(mod(sum(cds_exons(:,2)-cds_exons(:,1)),3)==0) ;
        if ~isempty(genes(i).cds_exons{j}), if ~(isequal(cds_exons, genes(i).cds_exons{j}(:,1:2))), fprintf('orfs differ\n') ; end ; end ;
        genes(i).tis{j} = tis_pos ;
        genes(i).cds_exons{j} = cds_exons ;
        genes(i).utr5_exons{j} = utr5_exons ;
        genes(i).utr3_exons{j} = utr3_exons ;
        % make sure tss and cleave exist
        if genes(i).strand=='+'
          genes(i).cdsStop{j} = cdsStop_pos-3 ;
          genes(i).tss{j} = genes(i).exons{j}(1,1);
          genes(i).cleave{j} = genes(i).exons{j}(end,2) ;
          assert(genes(i).tis{j}==genes(i).cds_exons{j}(1,1)) ;
          assert(genes(i).cdsStop{j}==genes(i).cds_exons{j}(end,2)-3) ;
        else
          genes(i).cdsStop{j} = cdsStop_pos+3 ;
          genes(i).tss{j} = genes(i).exons{j}(end,2) ;
          genes(i).cleave{j} = genes(i).exons{j}(1,1) ;
          assert(genes(i).tis{j}==genes(i).cds_exons{j}(end,2))
          assert(genes(i).cdsStop{j}==genes(i).cds_exons{j}(1,1)+3);
        end ;
        
        % sanity check
        if genes(i).strand=='+'
          cds_exons(:,2) = cds_exons(:,2)-1;
        else
          cds_exons(:,1) = cds_exons(:,1)+1;
        end
        [genestr] = load_genomic(genes(i).chr, genes(i).strand, cds_exons(:,1), cds_exons(:,2), genome_info, 1);
        assert(isequal(genestr(1:3), 'atg')) ;
        assert(any(ismember({'taa', 'tag', 'tga'}, genestr(end-2:end)))) ;
        
        num_tis_shift_max = num_tis_shift_max + 1 ;
        num_stop_shift_max = num_stop_shift_max + 1 ;
        
        if verb>=2,
          fprintf('=> successfully determined open reading frame of length %i\n', cds_len) ;
        end ;
      else
        if verb>=2,
          fprintf('=> Problems mapping TIS and/or cdsStop\n') ;
        end ;
      end ;
    else
      if verb>=2,
        if isempty(cand) || isempty(idx),
          fprintf('=> Cannot find appropriate ORF\n') ;
        else
          fprintf('=> Cannot find long enough ORF (%i)\n', cds_len) ;
        end ;
      end ;
    end ;
  end
end ; %%loop over genes
if verb==1,
  fprintf('Done.\n') ;
end ;

if verb>=1,
  fprintf('\t\t=> successfully determined maximal ORF\t%i\n', num_tis_shift_max) ;
  fprintf('\t\t=> successfully determined maximal ORF\t%i\n', num_stop_shift_max) ;
  fprintf('\n\n')
end ;

%keyboard
