function splicegraph = convert_ss2exon(sitegraph,strand)
% function splicegraph = convert_ss2exon(sitegraph,strand)
%
% Convert a splicegraph from splice site based to exon based.
% see convert_exon2ss for more details.


num_nodes = size(sitegraph{2},1);
num_exons = num_nodes;

splicegraph{1} = zeros(2,num_exons);
sitetype = zeros(2,num_exons);
% construct the exons, keeping track of the locations
idx_col = 0;
for ix = 1:num_nodes
  for iy = ix+1:num_nodes
    if sitegraph{2}(ix,iy) == 3
      % arrange so that iy is right of ix
      assert(sitegraph{1}(2,ix) ~= sitegraph{1}(2,iy))
      if sitegraph{1}(2,ix) < sitegraph{1}(2,iy)
	left_ss = sitegraph{1}(2,ix);
	right_ss = sitegraph{1}(2,iy);
	left_type = sitegraph{1}(1,ix);
	right_type = sitegraph{1}(1,iy);
      else
	left_ss = sitegraph{1}(2,iy);
	right_ss = sitegraph{1}(2,ix);
	left_type = sitegraph{1}(1,iy);
	right_type = sitegraph{1}(1,ix);
      end
      idx_col = idx_col + 1;
      splicegraph{1}(1,idx_col) = left_ss;
      splicegraph{1}(2,idx_col) = right_ss;
      sitetype(1,idx_col) = left_type;
      sitetype(2,idx_col) = right_type;
    end
  end
end
num_exons = idx_col;
splicegraph{1} = splicegraph{1}(:,1:num_exons);
sitetype = sitetype(:,1:num_exons);


% construct the introns
splicegraph{2} = zeros(num_exons);
for ix = 1:num_nodes
  for iy = ix+1:num_nodes
    if sitegraph{2}(ix,iy) == 7
      [ix_intron_left,ix_intron_right] = find_match_splicesite(splicegraph,...
						  sitegraph,sitetype,strand,ix,iy);
      splicegraph{2}(ix_intron_left,ix_intron_right) = 1;
      splicegraph{2}(ix_intron_right,ix_intron_left) = 1;
    end
  end
end


function [ix_intron_left,ix_intron_right] = find_match_splicesite(splicegraph,...
						  sitegraph,sitetype,strand,ix,iy)
% type must be the same
%ix_intron_start = floor(find(splicegraph{1}==sitegraph{1}(2,ix))/2);
%ix_intron_end = ceil(find(splicegraph{1}==sitegraph{1}(2,iy))/2);

if strand == '+'
  ix_intron_left = find(splicegraph{1}==sitegraph{1}(2,ix));
  ix_intron_left(find(sitetype(ix_intron_left)~=3))=[];
  
  ix_intron_right = find(splicegraph{1}==sitegraph{1}(2,iy));
  ix_intron_right(find(sitetype(ix_intron_right)~=6))=[];
else
  ix_intron_left = find(splicegraph{1}==sitegraph{1}(2,ix));
  ix_intron_left(find(sitetype(ix_intron_left)~=6))=[];
  
  ix_intron_right = find(splicegraph{1}==sitegraph{1}(2,iy));
  ix_intron_right(find(sitetype(ix_intron_right)~=3))=[];
end
ix_intron_left = ix_intron_left/2;
ix_intron_right = (ix_intron_right+1)/2;


assert(~isempty(ix_intron_left) && ~isempty(ix_intron_right));
