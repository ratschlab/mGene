function model = define_addContentFeatures(model)
  


for nn = 1:length(model.track_names)
  name = sprintf('track_%i', nn);
  start = model.cnt_plifs+1 ;
  cnt=0;
  
  model.(name).intergenic = start+cnt ;cnt=cnt+1;
  model.(name).utr_exon   = start+cnt ;cnt=cnt+1;
  model.(name).cds_exon   = start+cnt ;cnt=cnt+1;
  model.(name).intron     = start+cnt ;cnt=cnt+1;
 
  if isfield(model.use.contents, 'rna_seq_polya')&&model.use.contents.rna_seq_polya
    model.(name).rna_seq_polya     = start+cnt ;cnt=cnt+1;
  end
  
  %cnt_name = ['cnt_' name]; 
  model.cnt_tracks(nn) = length(fieldnames(model.(name)));
  plif_ids_name = [name '_plif_ids'];
  model.(plif_ids_name) = [start:start+model.cnt_tracks(nn)-1];
  
  %% declare range
  %range_name = [name '_range'];
  %model.(range_name) = struct;
  %field_names = fieldnames(model.(name)) ;
  %for s = 1:length(field_names)
  %  model.(range_name) = setfield(model.(range_name), field_names{s},feature_ranges{nn});
  %end
  fields_name = [name '_fields'];
  model.(fields_name) = fieldnames(model.(name));
  model.cnt_plifs = model.cnt_plifs + model.cnt_tracks(nn);
end
