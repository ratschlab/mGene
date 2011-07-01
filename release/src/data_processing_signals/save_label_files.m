function save_label_files(trivial_regions, P,genes)
  
%  save_label_files(trivial_regions, PAR,genes)
%
%  for each regions finds all occuring consensus sites for a signal
%  defined in PAR.Signal_name
%  Determines label and coverage by ESTs, cDNA, anno. Saves to binary
%  files which can be read by interval_query.
%
% INPUT trivial_regions  
%             (if several regions have same chr_num and strand error! )
%       genes (should be merged genes, with transcripts mapped to annotation)
%  
%       PAR.FN.input_sig.fn_candsites
%       PAR.Signal_name
%       PAR.Signal.(Signal_name).label_fct
%       PAR.Signal.(Signal_name).consensus
%       PAR.Signal.(Signal_name).offset  
%       PAR.Signal.(Signal_name).lwin_big
%       PAR.Signal.(Signal_name).rwin_big
%       PAR.Signal.(Signal_name).type
%       PAR.Signal.(Signal_name).name
%       PAR.Signal.(Signal_name).Conf_names
%
%  see retrieve_signals, interval_query 

fn_candsites = P.fn_candsites;
label_source = P.label_source;
info = P.info;
if isfield( P,'Signal')
  Signal = P.Signal;
  type = 'signal';
elseif isfield( P,'Content')
  type = 'content';
  Signal = P.Content;
end

if nargin>3
  CHR = [genes.chr_num]; 
  STR = [genes.strand]; 
end

if length([trivial_regions.chr_num])>2*length(unique([trivial_regions.chr_num]))
  fprintf(1, '[chr_num, strand] in regions are not unique!\n')
  error('use trivial_regions')
end

for r = 1:length(trivial_regions)
  strand =  trivial_regions(r).strand;
  fprintf(1,'region: %i\r',trivial_regions.id)
  filename = sprintf('%scontig_%i%s',fn_candsites,trivial_regions(r).chr_num,strand);
  if fexist([filename '.pos']) 
    fprintf(1,'pos file already exists.\r')
    continue
    fprintf(1,'\n')
  end
  if label_source.from_fn_candsites
    if ~fexist([filename '.mat'])
      error('fn_can with labels is missing')
    end
    load(filename,'pos','label')
    Conf = [];
    Signal.Conf_names = {};
  else
    seq = trivial_regions(r).seq;
    assert(length(seq)== trivial_regions(r).stop- trivial_regions(r).start+1)
    %% GET ALL LABELS
    
    geneidx = find(CHR==trivial_regions(r).chr_num & STR ==strand);
    if ~isempty(geneidx)
      fprintf(1,'number of genes in region: %i\n',length(geneidx))
    else
      fprintf(1,'number of genes in region: 0\n')
    end
    [pos,label,Conf] = feval(Signal.label_fct,genes(geneidx),seq,Signal,info,strand);
    if 0
      reg = retrieve_signals(trivial_regions(r),fn_old,'acc',{'label','est_cov','cdna_cov'},which_pos);
      length(setdiff(reg.Signals.acc.contig_pos,pos))
      [tt,ii1,ii2] = intersect(reg.Signals.acc.contig_pos,pos);
      assert(all(label(ii2)==reg.Signals.acc.label(ii1)')) 
      % assert(all(Conf(4,ii2)==reg.Signals.acc.est_cov(ii1)'))   
      % assert(all(Conf(5,ii2)==reg.Signals.acc.cdna_cov(ii1)'))   
      % keyboard
    end
    save(filename, '-V7', 'pos','label')
  end
  
  if isequal(filename(1),'~')
	home=getenv('HOME');
	home = deblank(home);
	filename=sprintf('%s/%s',home,filename(2:end)); 
  end
  
  if isequal(type,'signal')
    save_score_pos(pos,[label; Conf], filename ,{'label', Signal.Conf_names{:} });
  else
    save_score_pos(pos(:,1)', [pos(:,2)' ; label; Conf], filename ,{'pos2','label', Signal.Conf_names{:} });
  end
  [tmp,sha1sum] = unix(sprintf('sha1sum %s.pos',filename));
  save([filename '_sha1sum'], 'sha1sum')
end %%% loop over regions 
save([fn_candsites 'PAR.mat'],'P')
fprintf(1, 'save_label_files: DONE');
