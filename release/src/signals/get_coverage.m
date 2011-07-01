function ALL = get_coverage(ALL,gene,true_splice,label_2,info,info_names)
% [COV1,COV2] = get_coverage(ALL,gene,true_splice,info_names)
% 
% 
% get sum over all coverage for each true_splicesite and certain transcript()
% 
% INPUT  - gene:  gene structure of length 1;  
%          required fields:transcript_info,transcript_valid,exons,exons_confirmed, is_alt
%        - true_splice_site:    vector of all true splice sites (either
%          donor or acceptor, not both) that occur in the
%          splicegraph of the gene.
%        - info_names: cell which tells, you which transcripts to look
%          at ('anno','cDNA',...)
%   
% OUTPUT - cov1  number of coverage by given transcript type (according to info_names) for each true splice site. 
%        - cov2  number of times skipped by given transcript type (according to info_names) for each true splice site. 

  
COV1 = ALL.COV1 ;
COV2 = ALL.COV2;
POS = ALL.POS;
  
genic_idx = find(POS>gene.start&POS<gene.stop);
genic_pos = POS(genic_idx); 
assert(all(all(COV2(:,genic_idx)==0)))
assert(all(all(COV1(:,genic_idx)==0)))

cov1 = zeros(size(info_names,1),length(true_splice));
cov2 = zeros(size(info_names,1),length(true_splice));

% idx_const = find(sum(label_2==1)==0);
const_exons = gene.splicegraph{1}(:,find(~gene.exon_in_alt_region));
[tmp,idx_const,tmp] = intersect(true_splice,const_exons(:));
assert(all(sum(label_2(:,idx_const)==1)==0))
for type = 1:size(info_names,1) 
  if isequal(info_names{type,1},'anno')
    cov1(type,:) = -inf;
    cov2(type,:) = -inf;
  end
  t_idx = [];
  for tt1 = 1:size(info_names,2)
    if ~isempty(info_names{type,tt1})
      t_idx = [t_idx find(gene.transcript_info==info.(info_names{type,tt1}))];  
    end
  end
  if isempty(t_idx)
    continue
  end
  for t = t_idx
    if ~gene.transcript_valid(t)
      continue
    end
    [cov1_t, cov2_t] = get_coverage4transcript(gene,t,true_splice,info);
    if gene.transcript_info(t)==info.anno
      cov1(type,:) = max(cov1(type,:),cov1_t);
      cov2(type,:) = max(cov2(type,:),cov2_t); 
      if ~isempty(idx_const)&&~all(isinf(cov2_t(idx_const)))
        warning('here seems to be bug that should be fixed ... ') ;
        if keyboard_allowed()
          keyboard
        end
      end
    else
      cov1(type,:) = cov1(type,:) +cov1_t;
      cov2(type,:) = cov2(type,:) +cov2_t;
      if ~isempty(idx_const)&&~all(cov2_t(idx_const)==0)
        if keyboard_allowed()
          keyboard
        end
      end
    end
  end %%% loop over transcripts
end %%% loop over transcript types
idx_anno = find(strcmp(info_names,'anno')) ;
not_anno = setdiff(idx_anno,[1:size(info_names,1)]);

ii_alt2 = find(any(label_2==1));
ii_ir2 = find(label_2(2,:)==1);
if ~isempty(ii_alt2) && ~all(sum(cov2(not_anno,ii_alt2))>0|any(~isinf(cov2(idx_anno,ii_alt2))))
  warning('no cov2')
  % hold off;
  % viewsplicegraph(gene);
  if ~isempty(ii_ir2)
    fprintf('intron retention')
    % for a=1:length(true_alt_ir)
    %    hold on ; plot([true_alt_ir(a) true_alt_ir(a)],[1 10],'b')
    % end
  else
    %keyboard
  end
  gene.id
  warning('keyboard')
  %keyboard
end


[ii,ii1,ii2] = intersect(POS,true_splice);
assert(all(ismember(ii1,genic_idx)))

COV1(:,ii1) = cov1(:,ii2);
COV2(:,ii1) = cov2(:,ii2);

%%% SET COV2 OF NEGATIVE GENIC POSITIONS to min over coverage of
% neigbouring truesplicesites

true_splice = true_splice(ii2);
cov1 = cov1(:,ii2);


genic_neg_pos = setdiff(genic_pos,true_splice);
genic_cov2  =   zeros(size(info_names,1),length(genic_neg_pos));
[temp,ii1,ii2] = intersect(POS,genic_neg_pos);
assert(length(temp)==length(genic_neg_pos));
genic_neg_pos = genic_neg_pos(ii2);
if ~isequal(determine_engine, 'octave'),
  assert(issorted(true_splice))
end ;
for a=2:length(true_splice)
  idx = find(genic_neg_pos<true_splice(a)&genic_neg_pos>true_splice(a-1));
  genic_cov2(:,idx) = repmat(min([cov1(:,a),cov1(:,a-1)]'),length(idx),1)';
end
assert(all(all(COV2(:,ii1)==0)))
assert(all(all(COV1(:,ii1)==0)))
COV2(:,ii1) =genic_cov2;
    






function [cov1, cov2] = get_coverage4transcript(gene,transcript_idx,true_splice,info)

% [cov1, cov2] = get_coverage(gene,transcrit_idx,true_splice)
%
% for a given transcript transcrit_idx of the gene counts number of coverage for each
% of the true splice sites. 
% returned is cov1 the number of times a splice site is used in one
% evidence (e.g. transcript). In the case of ESTs this can be more than
% one, because a single transcript is produced out of a cluster of
% several ESTs.This information is stored in the
% field exons_confirmed{transcript_idx}, which is a Nx2 matrix, with the
% first number giving the confirmation for the exon and the second number
% giving the confirmation for the succeding
% intron. (exons_confirmed{transcript_idx}(end,2) shoeld always be 0, as
% there is no intron after the last exon.)
% The second returned value cov2 gives the number of evidences that
% disagree with the splicesite. e.i. the number of times a splice site is skipped.
   
  

exons = gene.exons{transcript_idx}(:,1:2);
exons_conf = gene.exons_confirmed{transcript_idx};

cov1 = -inf(1,length(true_splice));
cov2 = -inf(1,length(true_splice));

for a = 1:length(true_splice)
  %%% CHECK COVERAGE OF INCLUDED FORM     
  [idx1,idx2] = find(exons==true_splice(a));
  if ~isempty(idx1)
    assert(length(idx1)==1)
    cov1(a) = exons_conf(idx1,1); 
  end
  %%% CHECK EXON COVERAGE OF SKIPPED FORM    
  idx  = find(exons(:,1)<true_splice(a)&exons(:,2)>true_splice(a)) ;
  if ~isempty(idx)   
    assert(length(idx)==1)
    % cov2(a) = exons_conf(idx,2); 
    if gene.transcript_info(transcript_idx)==info.est|gene.transcript_info(transcript_idx)==info.cDNA &&...
          (idx==1|idx==size(exons,1))
      if all(label_2(:,a)==-1)
        continue
      else
        % alternative 3' or 5' or first or last exon; est ends
        % are usually unconfirmed 
      end
    elseif gene.transcript_info(transcript_idx)==info.anno
      if min(abs(true_splice(a)-exons(idx,:)))<3 
        % keyboard
        continue
      end
    end
    cov2(a) = exons_conf(idx,1);
  end  
  %%% CHECK INTRON COVERAGE OF SKIPPED FORM
  idx  = find(exons(1:end-1,2)<true_splice(a)&exons(2:end,1)>true_splice(a)) ;
  if ~isempty(idx)            
    assert(length(idx)==1)           
    if gene.transcript_info(transcript_idx)==info.anno
      cov2(a) = max(cov2(a),exons_conf(idx,2));      
    else
      cov2(a) = cov2(a)+exons_conf(idx,2);
    end
  end 
end %% loop over true_splice sites

assert(~any(~isinf(cov1)&~isinf(cov2)))




