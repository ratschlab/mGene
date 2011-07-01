function blocks_return = make_chunks(blocks)

%chunksize = 1.5e7;
chunksize = 5e6;
overlap = 5e5;
%warning('bad chunk size')
%chunksize = 1e6;
%overlap = 2e5;
blocks_return = [];

for j = 1:length(blocks)
  if blocks(j).stop>2*chunksize
    chunks = chunksize:chunksize:blocks(j).stop;
    if blocks(j).stop-chunks(end)<chunksize/2
      chunks(end)=[];
    end
    new_blocks = blocks(j);
    new_blocks.stop = chunks(1)+overlap;
    for k = 1:length(chunks)-1
      new_blocks(k+1) = blocks(j);
      new_blocks(k+1).start = chunks(k)-overlap;
      new_blocks(k+1).stop = chunks(k+1)+overlap;
    end
    new_blocks(end+1) = blocks(j);
    new_blocks(end).start = chunks(end)-overlap;
  else 
    new_blocks = blocks(j);
  end
  blocks_return = [blocks_return new_blocks];
end

for j = 1:length(blocks_return)
    blocks_return(j).id = j;
    if blocks_return(j).strand=='+'
      blocks_return(j).offset = blocks_return(j).start - 1;
    else
      blocks_return(j).offset = blocks_return(j).stop + 1;
    end
end

