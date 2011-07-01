function ok_idx = check_block_model_consistency(blocks,model)
%function ok_idx = check_block_model_consistency(blocks,model)
% walk through block.features along the 'true' path and 
% check consistency with model

ok = zeros(1,length(blocks));
for b=1:length(blocks)
  if mod(b,10)==0, fprintf('%i/%i ... ', b, length(blocks)) ; end 
  ok(b) = check_block(blocks(b),model,b);
end
ok_idx = find(ok);
fprintf('Done.\n') ;

return

function ok = check_block(block,model,num)
seq = lower(block.seq);
ok=1;
for t=1:length(block.truth)
  for j=1:length(block.truth(t).path)-1
    from_state   = block.truth(t).path(j);
    to_state     = block.truth(t).path(j+1);
    from_sig_pos = block.truth(t).pos_idx(j);
    to_sig_pos   = block.truth(t).pos_idx(j+1);
    from_pos     = block.all_pos(block.truth(t).pos_idx(j));
    to_pos       = block.all_pos(block.truth(t).pos_idx(j+1));
    

    assert(from_pos<to_pos)
    assert(from_sig_pos<to_sig_pos)

    % check segment lengths
    if strcmp(model.states(from_state).name,'intergenic_long') && strcmp(model.states(to_state).name,'intergenic_long')
	%this special case is necessary because the transition is zero 
	lengthrange  = [0 inf] ;%model.lengths_range.intergenic_long;
    %elseif strcmp(model.states(to_state).name,'intergenic_long')
	%lengthrange  = model.lengths_range.intergenic;
      length_name = 'intergenic_long';
    else
      transition   = model.transition_pointers(to_state,from_state,1);%transiton_pointers is transposed compared to model.A
      fields = fieldnames(model.lengths);
      length_name = '' ;
      for j=1:length(fields), if model.lengths.(fields{j})==transition, length_name = fields{j}; break, end,end
      assert(~isempty(length_name)) ;
      lengthrange  = model.lengths_range.(length_name);
    end
    if (to_pos-from_pos<lengthrange(1) ||to_pos-from_pos>lengthrange(2))%&~isequal(length_name,'intergenic_long')
      fprintf(1,'\ndiscard block %i, because its length (%i) is out of bounds (%i,%i) for %s\n',num,to_pos-from_pos,lengthrange(1),lengthrange(2), length_name);
      ok=0;
      %keyboard
      return
    end
    % check don/acc_1b transitions: does splicing lead to in frame stop codons?
    if strcmp(model.states(from_state).name,'CDS_don_1b')
      assert(seq(from_pos-1)=='t');
      assert(strcmp(model.states(to_state).name,'CDS_acc_1b'));
      if ismember(seq(to_pos+[0,1]),model.states(to_state).non_consensus)
        fprintf(1,'\ndiscard block %i because of a in frame stop codon\n',num);
        ok=0;
        return
      end
    end
    % check don/acc_2b transitions: does splicing lead to in frame stop codons?
    if strcmp(model.states(from_state).name,'CDS_don_2b')
      assert(seq(from_pos-2)=='t');
      assert(seq(from_pos-1)=='a'||seq(from_pos-1)=='g');
      assert(strcmp(model.states(to_state).name,'CDS_acc_2b'));
      if ~ismember(seq(to_pos+[-2:0]),model.states(to_state).consensus)
        fprintf(1,'\ndiscard block %i because of a in frame stop codon\n',num);
        ok=0;
        return
      end
    end
    for con=1:length(model.states(to_state).consensus)
      eq(con) = strcmp(model.states(to_state).consensus{con},seq((to_pos:to_pos+length(model.states(to_state).consensus{con})-1)-model.states(to_state).consensus_pos(con)));
    end
    if length(model.states(to_state).consensus)>0 && ~any(eq)
      fprintf(1,'\ndiscard block %i because of wrong %s-consensus at %i\n',num,model.states(to_state).name ,to_pos);
      ok=0;
      return
    end
    if isinf(model.A(from_state,to_state))
      fprintf(1,'\ndiscard block %i because it uses forbidden transition from %i to %i\n',num,from_state,to_state);
      ok=0;
      % keyboard
      return
    end

    if isinf(block.features{1}(to_state,to_sig_pos))
      error('discard block %i: feature entry is inf; to_state: %i to_pos: %i\n',num,to_state,to_sig_pos);
    end
  end
end
return
