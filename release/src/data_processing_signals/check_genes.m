function [genes, valid_map] = check_genes(genes, genome_info, which_checks) ;
% genes = check_genes(genes, genome_info,which_checks) ;

% TO DO CHECK SPLICE GRAPH: 
%  CHECK if splice graoh is connected. SHOULd WORK NOW
%  
%addpath ~/svn/tools/genomes/

min_boundary_dist=20 ;

if nargin<3 | isempty(which_checks), 
  which_checks = struct;
end

if ~isfield(which_checks,'exons_sorted')
  which_checks.exons_sorted = 1;
end ;
if ~isfield(which_checks,'intron_length')
  which_checks.intron_length = 0;
end
if ~isfield(which_checks,'exon_length')
  which_checks.exon_length = 0;
end
if ~isfield(which_checks,'splicesites');
  which_checks.splicesites = 0;
end
if ~isfield(which_checks,'cmp_valid_map');
  which_checks.cmp_valid_map = [] ;
end
if ~isfield(which_checks,'orf')
  which_checks.orf = 1;
end
if ~isfield(which_checks,'gene_length')
   which_checks.gene_length = 0;
end
if ~isfield(which_checks,'graph')
  which_checks.graph = 0;
end ;
if ~isfield(which_checks,'transacc')
  which_checks.transacc = 0;
end ;
if ~isfield(which_checks,'complete')
  which_checks.complete = 1;
end ;
which_checks


num_invalid_transcripts = 0 ;
num_invalid_genes = 0 ;
num_no_exons = 0;
num_no_trans = 0;
num_splice_not_ok = 0 ;
num_intron_not_ok = 0 ;
num_orf_not_ok = 0 ;
num_tis_not_ok = 0 ;
num_cds_exon_not_ok = 0 ;
num_stop_not_ok = 0 ;
num_graphs_not_ok = 0 ;
num_transacc_not_ok = 0;
num_at_boundary = 0;
num_stop_not_consistent = 0;
num_tis_not_consistent = 0;
genes_not_ok = [];

valid_map=[] ;

fprintf('checking %i genes: ', length(genes)) ;
for i=1:length(genes), 
  if mod(i,1000)==0||i==length(genes) 
    %fprintf(['gene=%i  invalid_genes=%i  invalid_transcripts=%i  no_exons=%i no_transcripts=%i ' ...
    %         'num_splice_not_ok=%i cds_problems=%i,%i,%i,%i num_graphs_not_ok =%i num_trans_not_ok =%i \r'],...
    %        i,num_invalid_genes,num_invalid_transcripts , num_no_exons,num_no_trans,num_splice_not_ok,...
    %        num_tis_not_ok, num_cds_exon_not_ok, num_orf_not_ok,
    %        num_stop_not_ok, num_graphs_not_ok,
    %        num_transacc_not_ok); 
    fprintf('%i ... ', i) ;
  end
  transcript_valid = [] ;
  is_valid = 1 ;
  
  if isempty(genes(i).exons)||all(cellfun('isempty',genes(i).exons))
    genes(i).is_valid = 0;
    num_invalid_genes = num_invalid_genes+1;
    num_no_exons = num_no_exons+1;
    continue
  end
  if isempty(genes(i).transcripts)||all(cellfun('isempty',genes(i).transcripts))
    genes(i).is_valid = 0;
    num_no_trans = num_no_trans+1;
    num_invalid_genes = num_invalid_genes+1;
    continue
  end
  
  
  % CHECK TRANSACC CONSENSUS
  trans_not_ok = 0;
  rm_trans=[];
  if which_checks.transacc 
    [transacc,idx] = unique(genes(i).transacc);
    for trans_=transacc'
      [mRNA, tmp, ok] =  load_genomic(genes(i).chr, char(genes(i).strand),...
                                      trans_-5, trans_+5,genome_info,1);
      if ~ok,
        warning('load genomic failed') ;
      else
        if ~isequal(mRNA(4:5),'ag')% &genes(i).strand=='-';
          mRNA(2:6)
          idx=find(ismember(genes(i).transacc,trans_));
          genes(i).transacc_info
          trans_not_ok = trans_not_ok+1;
          rm_trans=[rm_trans trans_];
        end
      end ;
    end
    idx=find(ismember(genes(i).transacc,rm_trans));
    if trans_not_ok
      warning('trans_not_ok');
      %keyboard
    end
    genes(i).transacc(idx)=[];
    genes(i).transacc_info(idx)=[];
    genes(i).transacc_conf(idx)=[];
    num_transacc_not_ok = num_transacc_not_ok+trans_not_ok; 
    
  end

  for j=1:length(genes(i).exons),


    [tmp,idx] = sort(genes(i).exons{j}(:,1)) ;
    
    genes(i).exons{j} = genes(i).exons{j}(idx,:) ;
    %if isfield(genes, 'exons_confirmed') && ~isempty(genes(i).exons_confirmed) && length(genes(i).exons_confirmed)>=j && ~isempty(genes(i).exons_confirmed{j})
    %  genes(i).exons_confirmed{j} = genes(i).exons_confirmed{j}(idx,:) ;
    %end
    

    %%% CHECK BOUNDARY
    % keyboard
    boundary_not_ok=0;
    if isempty(genes(i).start) || isempty(genes(i).stop)
      boundary_not_ok=1;
    else
      exons = genes(i).exons{j};
      if any(any(exons(:,1:2)<genes(i).start-1|exons(:,1:2)>genes(i).stop+1))
        warning('overlaying transcript found');
        boundary_not_ok=1;
      end
    end ;
   
    %%% CHECK EXONS SORTED
    exons_not_ok = 0 ;
    if isempty(genes(i).exons{j}),
      exons_not_ok = 1 ;
    else
      if which_checks.exons_sorted
        [tmp,idx] = sort(genes(i).exons{j}(:,2)) ;
        if ~isequal(idx',1:length(idx)), 
          fprintf('exons not sorted\n') ;
          exons_not_ok = 1 ;
        end ;
      end
    end ;
    %%% CHECK EXON LENGTH
    if which_checks.exon_length
      exon_length = genes(i).exons{j}(:,2)-genes(i).exons{j}(:,1);
      if any(exon_length<10)
        exons_not_ok = 1;
      end ;
    end

    
    %%% SPLICE SITE CONSENSUS
    splice_not_ok = 0 ;  
    
    if which_checks.splicesites && size(genes(i).exons{j},1)>1 && ~isempty(genes(i).exons{j}),
      % if any(genes(i).exons{j}(:,1)>genes(i).exons{j}(:,2))
      %    valid_splice(1:2,1:2)=0
      % else
      exons=genes(i).exons{j} ;
      [mRNA, valid_splice, ok] =  load_genomic(genes(i).chr, char(genes(i).strand), exons(:,1), exons(:,2),genome_info,1);
      if isfield(which_checks.cmp_valid_map, 'valid_splice'),
          valid_splice = valid_splice | (~which_checks.cmp_valid_map(i,j).valid_splice(:,1:2)) ;
      end ;
      if nargout>=2,
          if ~isempty(genes(i).cds_exons{j}),
              valid_splice_cds = (genes(i).exons{j} >= genes(i).cds_exons{j}(1,1)) & (genes(i).exons{j} <= genes(i).cds_exons{j}(end,2)) ;
          else
              valid_splice_cds = zeros(size(genes(i).exons{j})) ;
          end ;
          valid_map(i,j).valid_splice = [valid_splice valid_splice_cds] ;
      end ;

      if ~ok,
          splice_not_ok = 1 ;
      else
          % end
          if genes(i).strand=='+'
              splice_not_ok = sum(sum(valid_splice(2:end-1,:)')<2) + ~valid_splice(1,2)+~valid_splice(end,1);
          else
              splice_not_ok = sum(sum(valid_splice(2:end-1,:)')<2) + ~valid_splice(1,1)+~valid_splice(end,2);
          end
      end ;
      if splice_not_ok
          genes_not_ok = [genes_not_ok i];
      end
      num_splice_not_ok = num_splice_not_ok+splice_not_ok; 
    end
    
    %%% CHECK INTRON LENGTH
    intron_not_ok = 0 ;
    if which_checks.intron_length
      intron_length = genes(i).exons{j}(2:end,1)-genes(i).exons{j}(1:end-1,2);
      if any(intron_length>genome_info.max_intron_len)|...      
            any(intron_length<genome_info.min_intron_len)
        
        % fprintf('\nintron lengths: %i - %i \n',min(intron_length),max(intron_length)) ;
        intron_not_ok = sum(intron_length>genome_info.max_intron_len |...
                            intron_length<genome_info.min_intron_len);
      end ;
      num_intron_not_ok = num_intron_not_ok+intron_not_ok; 
    end
    transcript_valid(j) = ~any([exons_not_ok splice_not_ok intron_not_ok boundary_not_ok]);
  end %%loop over exons

  %%% CHECK ORF
  for j=1:length(genes(i).cds_exons), 
      if isfield(genes(i), 'transcript_coding'),
          if genes(i).transcript_coding(j)==0, continue ; end ;
      end ;
    if which_checks.orf ,
      cds_exons = genes(i).cds_exons{j};
      if isempty(cds_exons), 
        transcript_valid(j)=0;
        continue; 
      end
      if all(cds_exons(:)==0),
         genes(i).cds_exons{j}=[] ;
         transcript_valid(j)=0;
         continue ;
      end ;
      if genes(i).strand=='+'
        if ~isempty(genes(i).tis)
          if ~(genes(i).tis{j}==cds_exons(1,1)),
			num_tis_not_consistent = num_tis_not_consistent+1;
            transcript_valid(j) =0;
            continue ;
          end ;
        end
        if ~isempty(genes(i).cdsStop)
          if ~(genes(i).cdsStop{j}==cds_exons(end,2)-3),
			num_stop_not_consistent = num_stop_not_consistent+1;
            transcript_valid(j) =0;
            continue ;
          end ;
        end
        cds_exons(:,2) = cds_exons(:,2)-1;
      else
        if ~isempty(genes(i).tis)
          if ~(genes(i).tis{j}==cds_exons(end,2))
			num_tis_not_consistent = num_tis_not_consistent+1;
            transcript_valid(j) =0;
            continue ;
          end
        end
        if ~isempty(genes(i).cdsStop) 
          if ~(genes(i).cdsStop{j}==cds_exons(1,1)+3)
			num_stop_not_consistent = num_stop_not_consistent+1;
            transcript_valid(j) =0;
            continue ;
          end
        end
        cds_exons(:,1) = cds_exons(:,1)+1;
      end

      [cds, valid_splice, ok] =  load_genomic(genes(i).chr, char(genes(i).strand), cds_exons(:,1), cds_exons(:,2), genome_info, 1);
      if ~ok,
        transcript_valid(j) =0;
        continue ;
      end ;
      f = find_open_frames(cds,1) ;

      % keyboard
      % if first (last) cds_exon is smaler than 3 the start (stop) codon is 
      % spliced together, thus we will not have predictions for it
      too_short = cds_exons(1,2)-cds_exons(1,1)+1<3 || cds_exons(end,2)-cds_exons(end,1)+1<3;

      if too_short,
          if ~isfield(which_checks.cmp_valid_map, 'valid_cds_exon_len') || which_checks.cmp_valid_map(i,j).valid_cds_exon_len==1,
              num_cds_exon_ok = num_cds_exon_not_ok + 1 ;
              transcript_valid(j) =0;
              if nargout>=2, valid_map(i,j).valid_cds_exon_len = 0 ; end ;
          else
              if nargout>=2, valid_map(i,j).valid_cds_exon_len = 1 ; end ;
          end ;
      else
          if nargout>=2, valid_map(i,j).valid_cds_exon_len = 1 ; end ;
      end ;
      
      if too_short || ~any(ismember({'taa', 'tag', 'tga'}, cds(end-2:end)))
          if ~isfield(which_checks.cmp_valid_map, 'valid_stop') || which_checks.cmp_valid_map(i,j).valid_stop==1,
              transcript_valid(j) =0;
              num_stop_not_ok = num_stop_not_ok + 1;
              if nargout>=2, valid_map(i,j).valid_stop = 0 ; end ;
          else
              if nargout>=2, valid_map(i,j).valid_stop = 1 ; end ;
          end ;
      else
          if nargout>=2, valid_map(i,j).valid_stop = 1 ; end ;
      end ;
      if too_short || any(cds(1:3)~='atg'),
          if ~isfield(which_checks.cmp_valid_map, 'valid_tis') || which_checks.cmp_valid_map(i,j).valid_tis==1,
              transcript_valid(j) =0;
              num_tis_not_ok = num_tis_not_ok + 1;
              if nargout>=2, valid_map(i,j).valid_tis = 0 ; end ;
          else
              if nargout>=2, valid_map(i,j).valid_tis = 1 ; end ;
          end ;
      else
          if nargout>=2, valid_map(i,j).valid_tis = 1 ; end ;
      end ;

      if mod(length(cds),3)~=0 || f(1)~=1 || isempty(cds) 
          if ~isfield(which_checks.cmp_valid_map, 'valid_orf') || which_checks.cmp_valid_map(i,j).valid_orf==1,
              num_orf_not_ok  = num_orf_not_ok + 1 ;
              transcript_valid(j) =0;
              if nargout>=2, 
                  valid_map(i,j).valid_orf = 0 ; 
                  valid_map(i,j).valid_orf_info.cds_len = length(cds) ;
                  valid_map(i,j).valid_orf_info.cds_stop_pos = sort([strfind(lower(cds),'taa') strfind(lower(cds),'tag') strfind(lower(cds),'tga')]) ;
              end ;
          else
              if nargout>=2, valid_map(i,j).valid_orf = 1 ; end ;
          end ;
      else
          if nargout>=2, valid_map(i,j).valid_orf = 1 ; end ;
      end ;
    end
  end ; %% loop over cds_exons

  num_invalid_transcripts = num_invalid_transcripts + sum(~transcript_valid);

  % if isfield(genes(i),'transcript_valid') && ~isempty(genes(i).transcript_valid)
  %   genes(i).transcript_valid = genes(i).transcript_valid.*transcript_valid ;  
  % else
  genes(i).transcript_valid = transcript_valid ;
  % end
    
  %%% check genes length (only the region of the valid transcripts)
  length_not_ok = 0;
  if which_checks.gene_length,
    starts = [] ;
    ends = [] ;
    for j=1:length(genes(i).exons),
      if genes(i).transcript_valid(j)&&~isempty(genes(i).exons{j})
        starts(end+1)= genes(i).exons{j}(1,1) ;
        ends(end+1)  = genes(i).exons{j}(end,2) ;
      end ;
    end ;
    if max(ends)-min(starts)>genome_info.max_intron_len*10,
      % fprintf('gene length: %i\n', max(ends)-min(starts)) ;
      length_not_ok = 1;
    end ;
  end
  
  %%% check splicegraph
  %% please document me ... (some more) !!
  graph_not_ok = 0;
  if which_checks.graph && isfield(genes, 'splicegraph'), 
    edges = genes(i).splicegraph{2} ;
    if ~isequal(edges,0)
      ss = 0 ;
      for e=1:size(edges,1), 
        ss = ss + edges^e ; 
      end ; 
      if any(ss==0)
        graph_not_ok = 1;
        % viewsplicegraph(genes(i))
        % keyboard
      end
      num_graphs_not_ok = num_graphs_not_ok+graph_not_ok ;
    end
  end
  % if ~isfield(genes(i),'is_valid') || isempty(genes(i).is_valid)
  %   genes(i).is_valid = 1;
  % end    
  genes(i).is_valid = any(genes(i).transcript_valid)&&~length_not_ok&&...
      ~graph_not_ok&& ~isempty(genes(i).exons);
  % genes(i).is_valid = genes(i).is_valid&any(genes(i).transcript_valid)&~length_not_ok&...
  %     ~graph_not_ok& ~isempty(genes(i).exons);
  
  num_invalid_genes = num_invalid_genes + ~genes(i).is_valid; 
  
end ; %%loop over genes
fprintf('Done.\n') ;

genes(1).transcript_complete = [];
if which_checks.complete
  for i=1:length(genes),
    if ~any(genes(i).transcript_valid==1),
      genes(i).is_complete = 0 ;
      continue ;
    end ;
    for k=1:length(genes(i).transcripts)
      complete = 1;

% does not work since exons can be split into a cds and a utr part      
%      intersection = intersect(genes(i).cds_exons,genes(i).exons,'rows');
%      if ~(size(intersection,1)==size(genes(i).cds_exons,1)
%        error('adsfadfasdf')
%      end
%      intersection = intersect(genes(i).utr3_exons,genes(i).exons,'rows');
%      if ~(size(intersection,1)==size(genes(i).utr3_exons,1)
%        error('adsfadfasdf')
%      end
%      intersection = intersect(genes(i).utr5_exons,genes(i).exons,'rows');
%      if ~(size(intersection,1)==size(genes(i).utr5_exons,1)
%        error('adsfadfasdf')
%      end
%      clear intersection

      if isempty(genes(i).exons)|| isempty(genes(i).exons{k})
        genes(i).transcript_complete(k) = 0;
        continue
      end
      if length(genes(i).cds_exons)<k||isempty(genes(i).cds_exons{k})
        complete = 0;
      end
      if isfield(genes(i), 'utr_3prime'),
          if length(genes(i).utr_3prime)<k||isempty(genes(i).utr_3prime{k}) || isempty(genes(i).cds_exons{k})
              complete = 0;
          elseif genes(i).strand=='+'&&~all(all(genes(i).utr_3prime{k}(:, 1:2)>=max(max(genes(i).cds_exons{k}(:, 1:2)))))
              complete = 0;
          elseif genes(i).strand=='-'&&~all(all(genes(i).utr_3prime{k}(:, 1:2)<=min(min(genes(i).cds_exons{k}(:, 1:2)))))
              complete = 0;
          end
      end ;
      if isfield(genes(i), 'utr3_exons'),
          if length(genes(i).utr3_exons)<k||isempty(genes(i).utr3_exons{k}) || isempty(genes(i).cds_exons{k})
              complete = 0;
          elseif genes(i).strand=='+'&&~all(all(genes(i).utr3_exons{k}(:, 1:2)>=max(max(genes(i).cds_exons{k}(:, 1:2)))))
              complete = 0;
          elseif genes(i).strand=='-'&&~all(all(genes(i).utr3_exons{k}(:, 1:2)<=min(min(genes(i).cds_exons{k}(:, 1:2)))))
              complete = 0;
          end
      end ;
      if isfield(genes(i), 'utr_5prime'),
          if length(genes(i).utr_5prime)<k||isempty(genes(i).utr_5prime{k}) || isempty(genes(i).cds_exons{k})
              complete = 0;
          elseif genes(i).strand=='+'&&~all(all(genes(i).utr_5prime{k}(:, 1:2)<=min(min(genes(i).cds_exons{k}(:, 1:2)))))
              complete = 0;
          elseif genes(i).strand=='-'&&~all(all(genes(i).utr_5prime{k}(:, 1:2)>=max(max(genes(i).cds_exons{k}(:, 1:2)))))
              complete = 0;
          end
      end ;
      if isfield(genes(i), 'utr5_exons'),
          if length(genes(i).utr5_exons)<k||isempty(genes(i).utr5_exons{k}) || isempty(genes(i).cds_exons{k})
              complete = 0;
          elseif genes(i).strand=='+'&&~all(all(genes(i).utr5_exons{k}(:, 1:2)<=min(min(genes(i).cds_exons{k}(:, 1:2)))))
              complete = 0;
          elseif genes(i).strand=='-'&&~all(all(genes(i).utr5_exons{k}(:, 1:2)>=max(max(genes(i).cds_exons{k}(:, 1:2)))))
              complete = 0;
          end
      end ;
      if ~isempty(genes(i).tis) && length(genes(i).tis)>=k && ~isempty(genes(i).tis{k})
        if genes(i).tis{k}<genes(i).start||genes(i).tis{k}>genes(i).stop
          complete = 0;
        end
      else
        complete = 0;
      end

      if complete && ~isempty(genes(i).cdsStop)&&length(genes(i).cdsStop)>=k &&~isempty(genes(i).cdsStop{k})
        if genes(i).cdsStop{k}<genes(i).start||genes(i).cdsStop{k}>genes(i).stop
          complete = 0;
        end
      else
        complete = 0;
      end
    
      if complete && ~isempty(genes(i).tss)&& length(genes(i).tss)>=k&&~isempty(genes(i).tss{k})
        if genes(i).tss{k}==genes(i).tis{k}||genes(i).tss{k}<genes(i).start||genes(i).tss{k}>genes(i).stop
          complete = 0;
        end
      else 
        complete = 0;
      end
      
      if complete && ~isempty(genes(i).cleave)&& length(genes(i).cleave)>=k && ~isempty(genes(i).cleave{k})
        if genes(i).cleave{k}==genes(i).cdsStop{k} || genes(i).cleave{k}<genes(i).start || genes(i).cleave{k}>genes(i).stop
          complete = 0;
        end
      else 
        complete = 0;
      end
      
      for j=1:length(genes(i).exons),
        if ~(all(all(genes(i).exons{j}(:,1:2)>=genes(i).start& ...
                     genes(i).exons{j}(:,1:2)<=genes(i).stop+1)))
          genes(i).start = min(min(genes(i).exons{j}));
          genes(i).stop = max(max(genes(i).exons{j}));
        end
      end

      gene_contig_len=contig_len(genes(i).chr, genome_info) ;
      if genes(i).start<min_boundary_dist || genes(i).stop>gene_contig_len-min_boundary_dist,
        genes(i).is_valid=0 ;
        num_at_boundary = num_at_boundary + 1 ;
      end ;
      
      genes(i).transcript_complete(k) = complete & genes(i).transcript_valid(k);
    end
    if ~any(genes(i).transcript_complete)
      %keyboard
    end
    genes(i).is_complete = any(genes(i).transcript_complete);
  end
end

fprintf('\nProblems found:\n') ;
	    fprintf('\ttotal invalid genes\t\t%i\n',num_invalid_genes) ;
	    fprintf('\twithout transcripts\t\t%i\n', num_no_trans) ;
	    fprintf('\ttotal invalid transcripts\t%i\n', num_invalid_transcripts) ;
	    fprintf('\twithout exons\t\t\t%i\n', num_no_exons) ;
	    fprintf('\tsplice consensus missing\t%i\n', num_splice_not_ok) ;
	    fprintf('\ttis consensus missing\t\t%i\n', num_tis_not_ok) ;
	    fprintf('\tcds exon too short\t\t%i\n', num_cds_exon_not_ok) ;
	    fprintf('\tORF missing\t\t\t%i\n', num_orf_not_ok) ;
	    fprintf('\tstop consensus missing\t\t%i\n', num_stop_not_ok) ;
	    fprintf('\tcdsStop not consistent with cds_exons\t\t%i\n', num_stop_not_consistent) ;
	    fprintf('\ttis not consistent with cds_exons\t\t%i\n', num_tis_not_consistent) ;
	    fprintf('\tgraphs not ok\t\t\t%i\n', num_graphs_not_ok) ;
	    fprintf('\ttransacceptor not ok\t\t%i\n', num_transacc_not_ok) ;
	    fprintf('\tgenes at contig boundary\t%i', num_at_boundary) ;

fprintf('\n\n')
