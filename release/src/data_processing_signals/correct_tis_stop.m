function genes = correct_tis_stop(genes, genome_info, min_orf_len, min_orf_sep, verb) ;
% genes = correct_tis_stop(genes, genome_info, min_orf_len, min_orf_sep, verb) ;

if nargin<3,
  min_orf_len = 300 ;
end ;
if nargin<4,
  min_orf_sep = 0.7 ;
end ;
if nargin<5,
  verb=1 ;
  verb=2 ;
end ;

num_tis_not_ok = 0 ;
num_tis_shift_first = 0 ;
num_tis_shift_max = 0 ;
num_stop_not_ok = 0 ;
num_stop_shift_3 = 0 ;
num_stop_shift_max = 0 ;
num_tis_stop = 0 ;

for i=1:length(genes), 
  if verb==1 && mod(i,1000)==0,
    fprintf('%i/%i ... ', i, length(genes)) ;
  end ;
    

  if isempty(genes(i).exons)||all(cellfun('isempty',genes(i).exons))
    continue
  end
  if isempty(genes(i).transcripts)||all(cellfun('isempty',genes(i).transcripts))
    continue
  end
  
  for j=1:length(genes(i).exons),
    %%% CHECK ORF
    for j=1:length(genes(i).cds_exons), 
      cds_exons = genes(i).cds_exons{j};
      if isempty(cds_exons), 
        continue; 
      end
      if all(cds_exons(:)==0),
        genes(i).cds_exons{j}=[] ;
        continue ;
      end ;
      if genes(i).strand=='+'
        %if ~isempty(genes(i).tis)
        %  assert(genes(i).tis{j}==cds_exons(1,1)) ;
        %end ;
        %if ~isempty(genes(i).cdsStop)
        %  assert(genes(i).cdsStop{j}==cds_exons(end,2)-3) ;
        %end
        cds_exons(:,2) = cds_exons(:,2)-1;
      else
        %if ~isempty(genes(i).tis)
        %  assert(genes(i).tis{j}==cds_exons(end,2))
        %end
        %if ~isempty(genes(i).cdsStop) 
        %  assert(genes(i).cdsStop{j}==cds_exons(1,1)+3);
        %end
        cds_exons(:,1) = cds_exons(:,1)+1;
      end
      
      [cds, tmp, ok] =  load_genomic(genes(i).chr, char(genes(i).strand), cds_exons(:,1), cds_exons(:,2), genome_info, 1);
      if ~ok,
        continue ;
      end ;
      num_tis_stop = num_tis_stop + 1 ;
      
      % check if stop codon is not included in cds exon
      if ~any(ismember({'taa', 'tag', 'tga'}, cds(end-2:end)))
        num_stop_not_ok = num_stop_not_ok + 1 ;  
        if genes(i).strand=='+'
          cds_exons(end, 2) = cds_exons(end, 2) + 3 ;
        else
          cds_exons(1,1) = cds_exons(1,1) - 3 ;
        end
        [cds] =  load_genomic(genes(i).chr, char(genes(i).strand), cds_exons(:,1), cds_exons(:,2), genome_info, 1);
        
        % messy cds-stop correction
        if any(ismember({'taa', 'tag', 'tga'}, cds(end-2:end)))
          if genes(i).strand=='+'
            if ~isempty(genes(i).utr3_exons{j}),
              num_stop_shift_3 = num_stop_shift_3 + 1 ;
              genes(i).cds_exons{j}(end,2) = genes(i).cds_exons{j}(end,2)+3;
              genes(i).cdsStop{j} = genes(i).cdsStop{j}+3;
              genes(i).utr3_exons{j}(1,1) = genes(i).utr3_exons{j}(1,1) + 3 ;
              if genes(i).utr3_exons{j}(1,2)-genes(i).utr3_exons{j}(1,1)<1,
                ii=find(genes(i).exons{j}(:,2)==genes(i).utr3_exons{j}(1,2));
                if length(ii)~=1,
                  if verb>=2,
                    fprintf('When trying to correct CDS region: did not find correct exon (+).\n') ;
                  end ;
                  % undo action
                  genes(i).cds_exons{j}(end,2) = genes(i).cds_exons{j}(end,2)-3;
                  genes(i).cdsStop{j} = genes(i).cdsStop{j}-3;
                  genes(i).utr3_exons{j}(1,1) = genes(i).utr3_exons{j}(1,1) - 3 ;
                else
                  if genes(i).exons{j}(ii,2) < genes(i).utr3_exons{j}(1,1),
                    if verb>=2,
                      fprintf('Found (corrected) stop codon does not fit into annotated exon (%i-%i, %i-%i, +)\nDo not extend CDS %s by stop codon.\n', ...
                              genes(i).exons{j}(ii,1), genes(i).exons{j}(ii,2), genes(i).utr3_exons{j}(1,1), genes(i).utr3_exons{j}(1,2), ...
                              genes(i).transcripts{j}) ;
                    end ;
                    % undo action
                    genes(i).cds_exons{j}(end,2) = genes(i).cds_exons{j}(end,2)-3;
                    genes(i).cdsStop{j} = genes(i).cdsStop{j}-3;
                    genes(i).utr3_exons{j}(1,1) = genes(i).utr3_exons{j}(1,1) - 3 ;
                  else
                    genes(i).utr3_exons{j}(1,:)=[] ;
                  end
                end ;
              end ;
            end ;
          else
            if ~isempty(genes(i).utr3_exons{j}),
              num_stop_shift_3 = num_stop_shift_3 + 1 ;
              genes(i).cds_exons{j}(1,1) = genes(i).cds_exons{j}(1,1)-3;
              genes(i).cdsStop{j} = genes(i).cdsStop{j}-3;
              genes(i).utr3_exons{j}(end,2) = genes(i).utr3_exons{j}(end,2) - 3 ;
              if genes(i).utr3_exons{j}(end,2)-genes(i).utr3_exons{j}(end,1)<1,
                ii=find(genes(i).exons{j}(:,1)==genes(i).utr3_exons{j}(end,1));
                if length(ii)~=1,
                  if verb>=2,
                    fprintf('When trying to correct CDS region: did not find correct exon (-).\n') ;
                  end ;
                  % undo action
                  genes(i).cds_exons{j}(1,1) = genes(i).cds_exons{j}(1,1)-3;
                  genes(i).cdsStop{j} = genes(i).cdsStop{j}-3;
                  genes(i).utr3_exons{j}(end,2) = genes(i).utr3_exons{j}(end,2) - 3 ;
                end ;
                if genes(i).exons{j}(ii,1) > genes(i).utr3_exons{j}(end,2),
                  if verb>=2,
                    fprintf('Found (corrected) stop codon does not fit into annotated exon (%i-%i, %i-%i, -)\nDo not extend CDS %s by stop codon.\n', ...
                            genes(i).exons{j}(ii,1), genes(i).exons{j}(ii,2), genes(i).utr3_exons{j}(end,1), genes(i).utr3_exons{j}(end,2), ...
                            genes(i).transcripts{j}) ;
                  end ;
                  % undo action
                  genes(i).cds_exons{j}(1,1) = genes(i).cds_exons{j}(1,1)-3;
                  genes(i).cdsStop{j} = genes(i).cdsStop{j}-3;
                  genes(i).utr3_exons{j}(end,2) = genes(i).utr3_exons{j}(end,2) - 3 ;
                end
                genes(i).utr3_exons{j}(end,:)= [] ;
              end ;
            end ;
          end
        end
      end
      
      %have_stop1 = any(ismember({'taa', 'tag', 'tga'}, cds(end-2:end))) ;
      
      cds_exons = genes(i).cds_exons{j};
      if genes(i).strand=='+'
        cds_exons(:,2) = cds_exons(:,2)-1;
      else
        cds_exons(:,1) = cds_exons(:,1)+1;
      end
      [cds] =  load_genomic(genes(i).chr, char(genes(i).strand), cds_exons(:,1), cds_exons(:,2), genome_info, 1);
      
      have_tis = all(cds(1:3)=='atg') ;
      have_stop = any(ismember({'taa', 'tag', 'tga'}, cds(end-2:end))) ;
      
      %assert(have_stop==have_stop1) ;
      
      if have_tis && have_stop,
        continue ;
      end ;

      if verb>=2,
        fprintf('CDS region of gene %s ', genes(i).name) ;
      end ;
      if ~have_tis,
        num_tis_not_ok = num_tis_not_ok + 1 ;  
        if verb>=2,
          fprintf('does not start with ATG') ;
        end ;
      end ;
      if ~have_tis && ~have_stop,
        if verb>=2,
          fprintf(' and ') ;
        end ;
      end ;
      if ~have_stop,
        if verb>=2,
          fprintf('does not end with TAA, TGA, or TAG.\n') ;
        end ;
      else
        if verb>=2,
          fprintf('.\n') ;
        end ;
      end ;

      %if isequal(genes(i).name, 'gene:AGAP004308'), keyboard; end ;

      if ~have_tis && have_stop,
        idx=strfind(cds, 'atg') ;
        idx=idx(mod(idx,3)==1) ;
        if ~isempty(idx),
          atg_pos=idx(1) ;
          if length(cds)-atg_pos<min_orf_len,
            if verb>=2,
              fprintf('=> Shifting TIS to first ATG would lead to very short ORF (%i)\n', length(cds)-atg_pos) ;
            end ;
            continue ;
          end ;

          tis_pos = map_rna_pos(genes(i).cds_exons{j}, genes(i).strand, atg_pos) ;
          if genes(i).strand=='+',
            [cds_exons, utr5_exons, utr3_exons, ok] = split_exons_cds(genes(i).exons{j}, tis_pos, genes(i).cdsStop{j}+3, genes(i).strand) ;
          else
            [cds_exons, utr5_exons, utr3_exons, ok] = split_exons_cds(genes(i).exons{j}, tis_pos, genes(i).cdsStop{j}-3, genes(i).strand) ;
          end ;
          if ok,
            genes(i).tis{j} = tis_pos ;
            genes(i).cds_exons{j} = cds_exons ;
            genes(i).utr5_exons{j} = utr5_exons ;
            genes(i).utr3_exons{j} = utr3_exons ;
            % make sure tss and cleave exist
            if genes(i).strand=='+'
              genes(i).tss{j} = genes(i).exons{j}(1,1);
              genes(i).cleave{j} = genes(i).exons{j}(end,2) ;
              assert(genes(i).tis{j}==genes(i).cds_exons{j}(1,1)) ;
              assert(genes(i).cdsStop{j}==genes(i).cds_exons{j}(end,2)-3) ;
            else
              genes(i).tss{j} = genes(i).exons{j}(end,2) ;
              genes(i).cleave{j} = genes(i).exons{j}(1,1) ;
              assert(genes(i).tis{j}==genes(i).cds_exons{j}(end,2))
              assert(genes(i).cdsStop{j}==genes(i).cds_exons{j}(1,1)+3);
            end ;
            %if isequal(genes(i).name, 'gene:AGAP004714'), keyboard ; end ;
            
            % sanity check
            if genes(i).strand=='+'
              cds_exons(:,2) = cds_exons(:,2)-1;
            else
              cds_exons(:,1) = cds_exons(:,1)+1;
            end
            [genestr] = load_genomic(genes(i).chr, genes(i).strand, cds_exons(:,1), cds_exons(:,2), genome_info, 1);
            assert(isequal(genestr(1:3), 'atg')) ;
            assert(any(ismember({'taa', 'tag', 'tga'}, genestr(end-2:end)))) ;

            num_tis_shift_first = num_tis_shift_first + 1 ;
            if verb>=2,
              fprintf('=> Shifting TIS to first ATG in annotated CDS region\n') ;
            end ;
          else            
            if verb>=2,
              fprintf('=> Problems mapping TIS and/or cdsStop\n') ;
            end ;
          end ;
          continue ;
        else
          % it is probably not fixable...
          if verb>=2,
            fprintf('=> No ATG found in annotated CDS region\n') ;
          end ;
          continue ;
        end ;
      end ;
      %if have_tis && ~have_stop,
      %  fprintf('=> Cannot find stop codon\n') ;
      %  continue ;
      %end ;

      exons = genes(i).exons{j};
      if genes(i).strand=='+'
        exons(:,2) = exons(:,2)-1;
      else
        exons(:,1) = exons(:,1)+1;
      end
      [genestr] = load_genomic(genes(i).chr, genes(i).strand, exons(:,1), exons(:,2), genome_info, 1);
      [atg, e]=find_max_orfs(genestr, 'atg') ;
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
        tis_pos     = map_rna_pos(genes(i).exons{j}, genes(i).strand, atg(cand(idx))) ;
        cdsStop_pos = map_rna_pos(genes(i).exons{j}, genes(i).strand, e(cand(idx))+3+1) ;
        
        [cds_exons, utr5_exons, utr3_exons, ok] = split_exons_cds(genes(i).exons{j}, tis_pos, cdsStop_pos, genes(i).strand) ;
        if ok,
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
            assert(genes(i).tis{j}==genes(i).cds_exons{j}(end,2));
            assert(genes(i).cdsStop{j}==genes(i).cds_exons{j}(1,1)+3);
          end ;

          % sanity check
          if genes(i).strand=='+'
            cds_exons(:,2) = cds_exons(:,2)-1;
          else
            cds_exons(:,1) = cds_exons(:,1)+1;
          end
          [genestr] = load_genomic(genes(i).chr, genes(i).strand, cds_exons(:,1), cds_exons(:,2), genome_info, 1);
          %keyboard

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
  end ; %% loop over cds_exons
end ; %%loop over genes
if verb==1,
  fprintf('Done.\n') ;
end ;

if verb>=1,
  fprintf('\n\nChecked %i TIS and cdsStop positions.\nProblems found:\n', num_tis_stop) ;
  fprintf('\ttis consensus missing\t\t\t\t%i\n', num_tis_not_ok) ;
  fprintf('\t\t=> corrected to first ATG in frame\t%i\n', num_tis_shift_first) ;
  fprintf('\t\t=> successfully determined maximal ORF\t%i\n', num_tis_shift_max) ;
  fprintf('\tstop consensus missing\t\t\t\t%i\n', num_stop_not_ok) ;
  fprintf('\t\t=> shifted by stop codon by 3nt\t\t%i\n', num_stop_shift_3) ;
  fprintf('\t\t=> successfully determined maximal ORF\t%i\n', num_stop_shift_max) ;
  fprintf('\n\n')
end ;

%keyboard
