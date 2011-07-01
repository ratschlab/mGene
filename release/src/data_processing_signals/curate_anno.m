function  genes = curate_anno(genes, genome_info, log_fname)
  % genes = curate_anno(genes, genome_info, log_fname)

if exist('log_fname')
  fd=fopen(log_fname, 'a') ;
end

all_trans_comp = 0;  
some_trans_comp = 0; 
add_utr3 = zeros(1,length(genes));
add_utr5 = zeros(1,length(genes));
add_cds = zeros(1,length(genes));

add_5utr_length = 10;
add_3utr_length = 10;

for i=1:length(genes),
  gene = genes(i);
  gene_contig_len = contig_len(gene.chr, genome_info) ;

  if ~gene.is_valid
    continue
  end
  idx = find(gene.transcript_complete) ;
  if all(gene.transcript_complete)
    all_trans_comp = all_trans_comp + 1;
  elseif ~isempty(idx) 
    % remove incomplete transcripts
    some_trans_comp = some_trans_comp + 1;
    gene.exons = gene.exons(idx) ;
    gene.exons_confirmed = gene.exons_confirmed(idx);
    gene.cds_exons = gene.cds_exons(idx) ;
    gene.transcripts = gene.transcripts(idx) ;
    gene.transcript_valid = gene.transcript_valid(idx) ;
    gene.transcript_complete = gene.transcript_complete(idx) ;
    
    if ~isempty(gene.transcript_status)
      gene.transcript_status = gene.transcript_status(idx) ;
    end  
    gene.utr5_exons = gene.utr5_exons(idx) ;
    gene.utr3_exons = gene.utr3_exons(idx) ;
    if ~isempty(gene.tis),
      gene.tis = gene.tis(idx) ;
    end ;
    if ~isempty(gene.cdsStop),
      gene.cdsStop = gene.cdsStop(idx) ;
    end ;
    gene.tss = gene.tss(idx) ;
    gene.cleave = gene.cleave(idx) ;
    try,
      gene.polya = gene.polya(idx) ;
    catch
      % warning('polya elements missing') ;
    end ;
    if isempty(idx),
      gene.is_valid = 0 ;
    end ;
  else
    % fill in missing entries for all transcripts
    for k=1:length(gene.transcripts)
      
      % exons must be annotated
      if length(gene.exons)<k||isempty(gene.exons{k})
        assert(gene.is_valid==0)
        continue
      end
      
      % check if cds_exons are annotated
      if length(gene.cds_exons)>=k||~isempty(gene.cds_exons{k})
        
        % cds exons and exons have identical start  and no utrs
        if isequal(gene.cds_exons{k}(1,1), gene.exons{k}(1,1))
          if gene.strand =='+' 
            % no 5' utrs annotated 
            if length(gene.utr5_exons)<k || isempty(gene.utr5_exons{k})
              if gene.exons{k}(1,1)>add_5utr_length,
                add_utr5(i) = add_utr5(i)+1; 
                gene.utr5_exons{k} = [gene.exons{k}(1,1) - add_5utr_length, gene.exons{k}(1,1)];
                gene.exons{k}(1,1) = gene.exons{k}(1,1) - add_5utr_length; 
                assert(gene.tis{k} == gene.cds_exons{k}(1,1));
                gene.tss{k} = gene.exons{k}(1,1);
              else
                gene.utr5_exons{k} = [gene.exons{k}(1,1), gene.exons{k}(1,1)];
                gene.transcript_valid(k) = 0 ;
              end ;
            else
              warning('cds exons and exons identical, utrs annoted, check parser')
            end
          else %gene.strand =='-'
               % no 3' utrs annotated 
            if length(gene.utr3_exons)<k || isempty(gene.utr3_exons{k})
  	      if gene.exons{k}(1,1)>add_3utr_length,
                add_utr3(i) = add_utr3(i)+1; 
                gene.utr3_exons{k} = [gene.exons{k}(1,1) - add_3utr_length gene.exons{k}(1,1)];
                gene.exons{k}(1,1) = gene.exons{k}(1,1) - add_3utr_length; 
                assert(gene.cdsStop{k} == gene.cds_exons{k}(1,1)+3);
                gene.cleave{k} = gene.exons{k}(1,1);
	      else
		gene.utr3_exons{k} = [gene.exons{k}(1,1) gene.exons{k}(1,1)];
                gene.transcript_valid(k) = 0 ;
              end ;
            else % cds_exons and exons are identical, but annotated utrs
              warning('cds exons and exons identical, utrs annoted, check parser')
            end
          end
          
        end
        
        % cds exons and exons have identical end and no utrs
        if isequal(gene.cds_exons{k}(end,2), gene.exons{k}(end,2))
          if gene.strand =='+'
            % no 3'utrs annotated
            if length(gene.utr3_exons)<k||isempty(gene.utr3_exons{k})
              if gene.exons{k}(end,2) + add_3utr_length<gene_contig_len
                add_utr3(i) = add_utr3(i)+1; 
                gene.utr3_exons{k} = [gene.exons{k}(end,2) gene.exons{k}(end,2) + add_3utr_length]; 
                gene.exons{k}(end,2) = gene.exons{k}(end,2) + add_3utr_length; 
                assert(gene.cdsStop{k} == gene.cds_exons{k}(end,2)-3);
                gene.cleave{k} = gene.exons{k}(end,2);
              else
		gene.utr3_exons{k} = [gene.exons{k}(end,2) gene.exons{k}(end,2)];
                gene.transcript_valid(k) = 0 ;
              end ;
            else
              warning('cds exons and exons identical, utrs annoted, check parser')
            end
          else %gene.strand =='-'
               % no 5' utrs annotated 
            if length(gene.utr5_exons)<k || isempty(gene.utr5_exons{k})
              if gene.exons{k}(end,2) + add_5utr_length<gene_contig_len,
                add_utr5(i) = add_utr5(i)+1; 
                gene.utr5_exons{k} = [gene.exons{k}(end,2) gene.exons{k}(end,2) + add_5utr_length]; 
                gene.exons{k}(end,2) = gene.exons{k}(end,2) + add_5utr_length; 
                assert(gene.tis{k} == gene.cds_exons{k}(end,2));
                gene.tss{k} = gene.exons{k}(end,2);
              else
                gene.utr5_exons{k} = [gene.exons{k}(end,2), gene.exons{k}(end,2)];
                gene.transcript_valid(k) = 0 ;
              end ;
            else % cds_exons and exons are identical, but annotated utrs
              warning('cds exons and exons identical, utrs annoted, check parser')
            end
          end %sif strand          
        end
        
        if isempty(gene.tss) || isempty(gene.tss{k}),
          fprintf('warning: adding missing tss site\n') ;
          if gene.strand == '+',
            gene.tss{k} = gene.exons{k}(1,1) ;
          else
            gene.tss{k} = gene.exons{k}(end,2) ;
          end ;
        end ;
        if isempty(gene.cleave) || isempty(gene.cleave{k}),
          fprintf('warning: adding missing cleavage site\n') ;
          if gene.strand == '+',
            gene.cleave{k} = gene.exons{k}(end,2)
          else
            gene.cleave{k} = gene.exons{k}(1,1) ;
          end ;
        end ;
        
        if gene.transcript_valid(k),
          try
            %  cds_exons and exons are not identical, transcript should be complete 
            if gene.strand =='+'
              assert(gene.tis{k} == gene.cds_exons{k}(1,1));
              assert(gene.cdsStop{k} == gene.cds_exons{k}(end,2)-3);
              assert(gene.tss{k} == gene.exons{k}(1,1));
              assert(gene.cleave{k} == gene.exons{k}(end,2));
              assert(gene.utr5_exons{k}(1,1)==gene.tss{k});
              assert(gene.utr5_exons{k}(end,2)==gene.tis{k});
              assert(gene.utr3_exons{k}(1,1)==gene.cdsStop{k}+3);
              assert(gene.utr3_exons{k}(end,2)==gene.cleave{k});
            else
              assert(gene.tis{k} == gene.cds_exons{k}(end,2));
              assert(gene.cdsStop{k} == gene.cds_exons{k}(1,1)+3);
              assert(gene.tss{k} == gene.exons{k}(end,2));
              assert(gene.cleave{k} == gene.exons{k}(1,1));
              assert(gene.utr5_exons{k}(1,1)==gene.tis{k});
              assert(gene.utr5_exons{k}(end,2)==gene.tss{k});
              assert(gene.utr3_exons{k}(1,1)==gene.cleave{k});
              assert(gene.utr3_exons{k}(end,2)==gene.cdsStop{k}-3);
            end
          catch
            warning('one of the final checks failed, invalidating transcript %s', gene.name) ;
            gene.transcript_valid(k)=0 ;
          end;
        end ;
      elseif length(gene.cds_exons)<k||isempty(gene.cds_exons{k})
        warning('find orf')
        add_cds(i) = add_cds(i)+1;
      end
    end % loop over trancripts     
    if gene.strand=='+'
      gene.start = min([gene.tss{:}]);
      gene.stop = max([gene.cleave{:}]);
    else
      gene.start = min([gene.cleave{:}]);
      gene.stop = max([gene.tss{:}]);
    end
    genes(i) = gene;
  end % if all transcripts are incomplete
  
end %%loop over genes

take_map=ones(1,length(genes)) ;
num_trans_removed = 0 ;
for i=1:length(genes),
  idx = find(genes(i).transcript_valid) ;
  if isempty(idx), take_map(i)=0 ; end ;
  num_trans_removed = num_trans_removed + length(genes(i).transcript_valid)-length(idx) ;
  
  genes(i).exons = genes(i).exons(idx) ;
  genes(i).exons_confirmed = genes(i).exons_confirmed(idx);
  genes(i).cds_exons = genes(i).cds_exons(idx) ;
  genes(i).transcripts = genes(i).transcripts(idx) ;
  genes(i).transcript_valid = genes(i).transcript_valid(idx) ;
  genes(i).transcript_complete = genes(i).transcript_complete(idx) ;
end ;
genes=genes(take_map==1) ;
num_gene_removed = sum(take_map==0) ;

fprintf(1, ' %i genes with all transcripts complete.\n', sum(all_trans_comp)) ;
fprintf(1, ' %i genes with at least one transcript complete.\n', sum(some_trans_comp)) ;
fprintf(1, ' %i 5UTRs added to transcripts of %i genes.\n',sum(add_utr5), sum(add_utr5>0) ) ;
fprintf(1, ' %i 3UTRs added to transcripts of %i genes.\n',sum(add_utr3), sum(add_utr3>0) ) ;
fprintf(1, ' %i max ORF determined for transcripts of %i genes.\n', sum(add_cds), sum(add_cds>0)) ;
fprintf(1, ' %i transcripts and %i genes removed at contig boundary\n', num_trans_removed, num_gene_removed) ;

if exist('fd')
fprintf(fd, '\n')
fprintf(fd, ' %i genes with all transcripts complete.\n', sum(all_trans_comp)) ;
fprintf(fd, ' %i genes with at least one transcript complete.\n', sum(some_trans_comp)) ;
fprintf(fd, ' %i 5UTRs added to transcripts of %i genes.\n',sum(add_utr5), sum(add_utr5>0) ) ;
fprintf(fd, ' %i 3UTRs added to transcripts of %i genes.\n',sum(add_utr3), sum(add_utr3>0) ) ;
fprintf(fd, ' %i max ORF determined for transcripts of %i genes.\n', sum(add_cds), sum(add_cds>0)) ;
fprintf(fd, ' %i transcripts and %i genes removed at contig boundary\n', num_trans_removed, num_gene_removed) ;
fclose(fd);
end
