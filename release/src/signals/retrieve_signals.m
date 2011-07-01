function regions = retrieve_signals(regions,base_dir,signal_name,fields,which_pos,offset)

% regions = retrieve_signals(regions,base_dir,signal_name,fields,which_pos)
%
% loads all candidate positions and addional informations as specified in
% fields,(i.e.'label','anno_cover','est_cover','cdna_cover','Conf_cum' )
% from binary fiels and adds them to the field Signals.signal_name of the
% structure regions 
%
% see also save_label_files
if nargin<6
  offset = 0;
end
  
if nargin<5
  which_pos.local = 1;
  which_pos.contig = 1;
end

genome_info = init_genome(regions(1).config);

f_name = sprintf('%scontig_%i%s.pos', base_dir,regions(1).chr_num,regions(1).strand );
if ~fexist(f_name)&&~isempty(base_dir),
  fprintf('file %s does not exist\n', f_name);
end

if ~isempty(fields) && ~isempty(base_dir)
  f_idx = [];
  for f=1:length(fields),
    f_name = sprintf('%scontig_%i%s.%s', base_dir,regions(1).chr_num,regions(1).strand, fields{f});
    if fexist(f_name),
      f_idx = [f_idx f];
    else
      fprintf('file %s does not exist\n', f_name);
    end
  end
  if ~isempty(f_idx),
    fields = fields(f_idx);
  else
    error('no files available')
  end
end

if ~isfield(regions,'Signals')
  regions(1).Signals = [];
end


for r=1:length(regions)
  if isfield(regions(r),'num_fragments')&&~isempty(regions(r).num_fragments)
    num_fragments = regions(r).num_fragments;
  else
    num_fragments = 1;
  end
  prefix = 0;
  
  SIGNALS.pos = [];
  SIGNALS.contig_pos = [];
  for f=1:length(fields)
    SIGNALS.(fields{f})= [];
  end
  if ismember('pos2',fields)
    SIGNALS.contig_pos2 = [];
  end
  for p=1:num_fragments
    contig_length = dir(genome_info.flat_fnames{regions(r).chr_num(p)}) ;
    contig_length = contig_length.bytes;
    %if 1%r==5145
    %  fprintf(1, 'regions.stop: %i, regions.start: %ir: %i',regions(r).start,regions(r).stop, r)
    %end
    if regions(r).start(p)<1
      warning('regions.start<1; set to 1!')
      regions(r).start(p) = 1;
    end
    if regions(r).stop(p)>contig_length
      warning('regions.stop>contig length; set to contig_length!')
      regions(r).stop(p) = contig_length;
    end
    f_name = sprintf('%scontig_%i%s',base_dir, regions(r).chr_num(p), ...
                     regions(r).strand(p));
    if fexist([f_name '_sha1sum.mat'])
      [tmp,sha1sum] = unix(sprintf('sha1sum %s.pos',f_name));
      sha1sum = separate(sha1sum,' ');
      Old = load([f_name '_sha1sum'],'sha1sum');
      Old.sha1sum = separate(Old.sha1sum,' ');
      if ~isequal(Old.sha1sum{1},sha1sum{1})
        warning('sha1sum different')
      end
      %  else
      %    warning('no sha1sum check performed')
    end
    
    if isempty(fields)
      fields{1} ='pos';
    end
    if regions(r).start(p)+offset>regions(r).stop(p)-offset
      regions(r).Signals.(signal_name) = SIGNALS;
      continue
    end
	if isempty(base_dir)
		% this is used to determine the contribution of the different signals 
		% to the gene finding score or the over all performance
		contig_pos = (regions(r).start(p):regions(r).stop(p))';
		score = zeros(length(contig_pos), length(fields));	
	else
    	d=dir([f_name '.pos']) ;
    	if d(1).bytes==0,
    	  contig_pos = [] ;
    	  score = zeros(0, length(fields)) ;
    	else
    	  [contig_pos,score] = interval_query(f_name,fields,[regions(r).start(p)+offset;regions(r).stop(p)-offset]);
    	end ;
    	if ismember('pos2',fields)
    	  if ~isempty(contig_pos)
    	    idx_p = find(strcmp('pos2',fields));
    	    contig_pos2 = score(:,idx_p);
    	  else
    	    contig_pos2 =[];
    	  end
    	end
	end
    assert(isequal(contig_pos, sort(contig_pos))) ;
    %assert(issorted(contig_pos))
    if ~isempty(contig_pos)
      if ~(contig_pos(1)>=regions(r).start(p)&contig_pos(end)<=regions(r).stop(p))
        if keyboard_allowed()
          keyboard
        end;
      end
    end
    if which_pos.local
      if regions(r).strand(p)=='+'
        pos = contig_pos - regions(r).start(p) + 1+prefix;
        if exist('contig_pos2', 'var')
          pos2 = contig_pos2 - regions(r).start(p) + 1+prefix;
        end
        SIGNALS.pos = [SIGNALS.pos pos];
      else
        pos = regions(r).stop(p) - contig_pos + 1;
        [pos,idx] = sort(pos);
        score = score(idx,:);
        contig_pos = contig_pos(idx);
        if exist('contig_pos2', 'var')
          contig_pos2 = contig_pos2(idx);
          pos2 = regions(r).stop(p) - contig_pos2 + 1;
          SIGNALS.pos = [SIGNALS.pos pos2];
        else
          SIGNALS.pos = [SIGNALS.pos pos];
        end
      end
    end
    if which_pos.contig
      SIGNALS.contig_pos = [SIGNALS.contig_pos contig_pos];
      if exist('contig_pos2', 'var')
        SIGNALS.contig_pos2 = [SIGNALS.contig_pos2 contig_pos2];
        if regions(r).strand(p)=='+'
          SIGNALS.pos2 = [SIGNALS.pos2 pos2];
        else
          SIGNALS.pos2 = [SIGNALS.pos2 pos];
        end
      end
    end
    for f=1:length(fields)
      if strmatch(fields{f},'label') ix=f; end
      if ~isequal(fields{f},'pos')&~isequal(fields{f},'pos2')
        SIGNALS.(fields{f})= [SIGNALS.(fields{f}), score(:,f)];
      end
    end
    if exist('ix', 'var')
      u = unique(score(:,ix));
      if length(u)==1 assert(u==-1|u==1),
      elseif length(u)==2 assert(all(u==[-1;1])),
      elseif isempty(u),
      else error('unique(label) is not [-1;1]');
      end
    end
    regions(r).Signals.(signal_name) = SIGNALS;
    prefix = prefix +regions(r).stop(p)- regions(r).start(p);
  end %%loop over fragments
end %% loop over blocks
