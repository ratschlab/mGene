function model = define_transition_pointers(model)
  
A = model.A;
mRNA_start_states_out=model.mRNA_start_states_out;
%%% transition_pointers contain for each allowed transition pointers to lengths plifs
transition_pointers = zeros(size(A)) ;
transition_pointers(~isinf(A))=A(~isinf(A)) ;
transition_pointers = transition_pointers' ;

%%% add pointers to content PLIFs
transition_pointers(:,:,2) = zeros(size(transition_pointers,1), size(transition_pointers,2), 1) ;
for i=1:size(model.plif_links,1)
  [idx1,idx2] = find(transition_pointers(:,:,1)==model.plif_links(i,1));
  for j=1:length(idx1)
    transition_pointers(idx1(j),idx2(j),2) = model.plif_links(i,2);
  end ;
end


%%% add pointers to linear feature PLIFs
names = {};
for j = 1:length(model.track_names)
  names{end+1} = sprintf('track_%i', j);
end
for j = 1:length(model.segment_feature_names);
  names{end+1} = sprintf('segment_feature_%i',j);
  names{end+1} = sprintf('segment_score_%i',j);
end

for j=1:length(names)
  name = names{j};
  link_name = [name '_links'];
  dim3 = size(transition_pointers,3);
  transition_pointers(:,:,dim3+1) = zeros(size(transition_pointers,1), size(transition_pointers,2), 1) ;
  for i=1:size(model.(link_name),1)
    [idx1,idx2] = find(transition_pointers(:,:,1)==model.(link_name)(i,1));
    for j=1:length(idx1)
      transition_pointers(idx1(j),idx2(j),dim3+1) = model.(link_name)(i,2);
    end ;
  end
end

%%% add pointers to frame PLIFs
tmp = zeros(size(transition_pointers(:,:,1))) ;
tmp(mRNA_start_states_out, model.state_ids.cdsStop(1)) = model.contents.cds_frame0 ;% single cds exon  

if ~isfield(model.use.contents, 'pre_comp') || ~model.use.contents.pre_comp,
  tmp(mRNA_start_states_out, model.state_ids.don(2)) = model.contents.cds_frame0 ;% first cds exon 
  tmp(mRNA_start_states_out, model.state_ids.don(3:4)) = model.contents.cds_frame1 ;% first cds exon  
  tmp(mRNA_start_states_out, model.state_ids.don(5:7)) = model.contents.cds_frame2 ;% first cds exon 
  
  tmp(model.state_ids.acc(2:7), model.state_ids.don(2)) = model.contents.cds_frame0 ;% middle cds exon 
  tmp(model.state_ids.acc(2:7), model.state_ids.don(3:4)) = model.contents.cds_frame1 ;% middle cds exon 
  tmp(model.state_ids.acc(2:7), model.state_ids.don(5:7)) = model.contents.cds_frame2 ;% middle cds exon 
  tmp(model.state_ids.acc(2:7), model.state_ids.cdsStop(1)) = model.contents.cds_frame0 ;% last cds exon 
else
  tmp(mRNA_start_states_out, model.state_ids.don(2:7)) = model.contents.cds_frame0 ;% first cds exon 
  tmp(model.state_ids.acc(2),model.state_ids.don(2:7)) = model.contents.cds_frame0 ;% middle cds exon 
  tmp(model.state_ids.acc(3:4),model.state_ids.don(2:7)) = model.contents.cds_frame1 ;% middle cds exon 
  tmp(model.state_ids.acc(5:7),model.state_ids.don(2:7)) = model.contents.cds_frame2 ;% middle cds exon 
  tmp(model.state_ids.acc(2), model.state_ids.cdsStop(1)) = model.contents.cds_frame0 ;% last cds exon 
  tmp(model.state_ids.acc(3:4), model.state_ids.cdsStop(1)) = model.contents.cds_frame1 ;% last cds exon 
  tmp(model.state_ids.acc(5:7), model.state_ids.cdsStop(1)) = model.contents.cds_frame2 ;% last cds exon 
end ;


dim3 = size(transition_pointers,3);
transition_pointers(:,:,dim3+1) = tmp' ; % don't touch this---tmp has to be transposed!


a_trans = define_a_trans(A, transition_pointers, model.seg_links);


model.transition_pointers = transition_pointers ;
model.a_trans = a_trans;
