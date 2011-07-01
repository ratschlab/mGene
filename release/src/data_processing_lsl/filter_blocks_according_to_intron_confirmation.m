function idx = filter_blocks_according_to_intron_confirmation(blocks)

if ~isfield(blocks, 'segment_lists')
  idx = 1:length(blocks);
  fprintf('no segment lists found\n')
  return
end

all_confirmed = zeros(1,length(blocks));
confirmation  = zeros(1,length(blocks));
num_introns   = zeros(1,length(blocks));
splice_conf   = zeros(1,length(blocks));

for j = 1:length(blocks)
  confirmed_introns = blocks(j).segment_lists{1};
  if isempty(blocks(j).truth(1).segments)
    annotated_introns = [];
  else
    annotated_introns = blocks(j).truth(1).segments(blocks(j).truth(1).segments(:,3)==5, 1:2);
  end
  if size(intersect(confirmed_introns, annotated_introns, 'rows'), 1)==size(annotated_introns, 1)
    all_confirmed(j) = 1;
  end
  confirmation(j) = size(intersect(confirmed_introns, annotated_introns, 'rows'), 1);
  num_introns(j) = size(annotated_introns, 1);
  num_conf(j) = length(unique(confirmed_introns));

  splice_conf(j) = length(intersect(unique(confirmed_introns) , unique(annotated_introns)));
end
cutoff = 1.4;

idx = find(num_conf./splice_conf<cutoff);

plot(splice_conf, num_conf, '.')
hold on
plot(1:100, (1:100)*cutoff)

fprintf('filter criterion: keep those blocks where at least %f%% of the splice sites \n from spliced reads correspond to annotated splice sites\n', 100/cutoff)
fprintf('keep %f%% of the blocks according to filter criterion\n', 100*mean(num_conf./splice_conf<cutoff))
