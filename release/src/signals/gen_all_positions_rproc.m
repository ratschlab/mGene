function [dummy1, dummy2]=gen_all_positions_rproc(P)

% gen_all_positions_rproc(P)
% P has fields:  
%   - region;  
%   - blocks 
%   - Signal 
%   - fn_pos 
%   - num_splits 
%   
  

dummy1=[] ;
dummy2=[] ;
  
paths  

region = P.region;  
blocks = P.blocks;
Signal = P.Signal;
fn_pos = P.fn_pos;
num_splits =  P.num_splits;


fprintf('start loading sequence for region %i%s\r',region.chr_num,region.strand)
tt = load_sequence(region);
if ~isempty(Signal.consensus)
  pos = find_consensus(tt.seq,Signal,region.strand);
elseif region.strand=='+'
  pos = Signal.lwin_big+1:10:length(tt.seq)-Signal.rwin_big; 
elseif region.strand=='-'
  pos = Signal.rwin_big+1:10:length(tt.seq)-Signal.lwin_big; 
end    

clear tt;

%-------------
% Determine SVM
%-------------
svm_num = zeros(size(pos));

% faster if there are many blocks and few splits
% needs more memory
if ~isempty(blocks),
  splits=[blocks.split] ; 
else
  splits=[] ;
end ;
assert(length(splits)==length(blocks)) ;
for split = unique(splits),
  fprintf('split %i\r',split)
  
  split_idx = find(split==splits) ;

  split_pos = {} ; p=length(split_idx); 
  for b_idx=split_idx,
    assert(~isfield(blocks, 'num_fragments') || blocks(b_idx).num_fragments<=1) ;
    if region.strand==blocks(b_idx).strand(1) && region.chr_num==blocks(b_idx).chr_num(1),
      if isfield(blocks(b_idx),'reg_no_overlap') && ~isempty(blocks(b_idx).reg_no_overlap)
        split_pos{p} = blocks(b_idx).start+[blocks(b_idx).reg_no_overlap(1):blocks(b_idx).reg_no_overlap(2)];
      else
        split_pos{p} = blocks(b_idx).start:blocks(b_idx).stop ;
      end
    end ;
    p=p-1 ;
  end ;
  split_pos = [split_pos{:}] ;
  [tmp,idx]=intersect(pos, split_pos); 
  assert(all(svm_num(idx)==0 | svm_num(idx)==split))
  svm_num(idx) = split ;
end ;

if 0, % slow version
  for b_idx = 1:length(blocks)
    fprintf('block %i\r',b_idx)
    block = blocks(b_idx);
    if ~isfield(block, 'num_fragments'),
      block.num_fragments = 1 ;
    end ;
    if ~isfield(block, 'reg_no_overlap'),
      block.reg_no_overlap = [block.start block.stop] ;
    end ;
    for f_idx = 1:block.num_fragments
      if f_idx>1
        error('block.reg_no_overlap not implemented for more than one fragment')
      end
      if region.strand==block.strand(f_idx) && region.chr_num==block.chr_num(f_idx);
        % idx =
        % find(pos>=block.start(f_idx)+PAR.regions.offset&pos<=block.stop(f_idx)-PAR.regions.offset);
        
        idx = find(pos>=block.reg_no_overlap(1) & pos<=block.reg_no_overlap(2));
        assert(all(svm_num(idx)==0 | svm_num(idx)==block.split))
        svm_num(idx) = block.split;
        % keyboard
      end
    end    
  end
end ;

filename = sprintf('%scontig_%i%s',fn_pos,region.chr_num,region.strand)
if isequal(filename(1),'~')
  filename = sprintf('%s/%s', getenv('HOME'), filename(2:end));
end

if ~isempty(find(svm_num==0))

  idx = find(svm_num==0);
  % nums = repmat([1:num_splits],1,ceil(length(idx)/num_splits));
  % nums = nums(randperm(length(nums)));
  % make blocks instead (faster computing for contents)

  num_blocks = round(max(length(idx)/10000, num_splits));
  nums = zeros(1, num_blocks*ceil(length(idx)/num_blocks)) ;
  p=1 ;
  for i=1:num_blocks
    nums( p : p+ceil(length(idx)/num_blocks)-1 ) = repmat( mod(num_blocks,num_splits) + 1, 1, ceil(length(idx)/num_blocks)) ;
    p = p + ceil( length(idx)/num_blocks ) ;
  end
  assert(p-1==length(nums));
  svm_num(idx) = nums(1:length(idx));

end

if ~isempty(pos)
  save_score_pos(pos,svm_num, filename ,{'svm'});
else
  unix(sprintf('touch %s.pos',filename));
  unix(sprintf('touch %s.svm',filename));
  warning('no position found')
end
