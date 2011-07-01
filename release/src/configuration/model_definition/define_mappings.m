function model = define_mappings(model, feature_names)
  

%%mappings between length and contents
  
model.plif_links = [] ;
model.plif_links(end+1,:) = [model.lengths.intergenic      model.contents.intergenic] ;
model.plif_links(end+1,:) = [model.lengths.first_cds_exon  model.contents.cds_exon] ;
model.plif_links(end+1,:) = [model.lengths.last_cds_exon   model.contents.cds_exon] ;
model.plif_links(end+1,:) = [model.lengths.middle_cds_exon model.contents.cds_exon] ;
model.plif_links(end+1,:) = [model.lengths.single_cds_exon model.contents.cds_exon] ;
model.plif_links(end+1,:) = [model.lengths.intron          model.contents.intron] ;
if model.use.segments.transexon
  model.plif_links(end+1,:) = [model.lengths.transexon       model.contents.utr_exon] ;
end
model.plif_links(end+1,:) = [model.lengths.utr5exon        model.contents.utr_exon] ;
model.plif_links(end+1,:) = [model.lengths.utr3exon        model.contents.utr_exon] ;
if model.use.segments.polya_tail
  model.plif_links(end+1,:) = [model.lengths.polya_tail      model.contents.utr_exon] ;
end
if model.use.segments.intercistronic
  model.plif_links(end+1,:) = [model.lengths.intercistronic  model.contents.intercistronic] ;
  model.plif_links(end+1,:) = [model.lengths.intergenictrans  model.contents.intergenic] ;
end
model.plif_links(end+1,:) = [model.lengths.intergenic_long model.contents.intergenic] ;

if isfield(model.use.segments, 'rna_seq_polya')&&model.use.segments.rna_seq_polya
  model.plif_links(end+1,:) = [model.lengths.rna_seq_polya  model.contents.utr_exon] ;
end
if isfield(model.use, 'non_coding') && model.use.non_coding
	model.plif_links(end+1,:) = [model.lengths.first_nc_exon  model.contents.nc_exon] ;
	model.plif_links(end+1,:) = [model.lengths.middle_nc_exon  model.contents.nc_exon] ;
	model.plif_links(end+1,:) = [model.lengths.last_nc_exon  model.contents.nc_exon] ;
end

%%mappings between length and additional contents
for j=1:length(model.track_names)
  name = sprintf('track_%i', j);
  link_name = [name '_links'];
  model.(link_name) = [] ;   
  model.(link_name)(end+1,:) = [model.lengths.intergenic      model.(name).intergenic] ;
  model.(link_name)(end+1,:) = [model.lengths.first_cds_exon  model.(name).cds_exon] ;
  model.(link_name)(end+1,:) = [model.lengths.last_cds_exon   model.(name).cds_exon] ;
  model.(link_name)(end+1,:) = [model.lengths.middle_cds_exon model.(name).cds_exon] ;
  model.(link_name)(end+1,:) = [model.lengths.single_cds_exon model.(name).cds_exon] ;
  model.(link_name)(end+1,:) = [model.lengths.intron          model.(name).intron] ;
  if model.use.segments.transexon
    model.(link_name)(end+1,:) = [model.lengths.transexon       model.(name).transexon] ;
  end
  model.(link_name)(end+1,:) = [model.lengths.utr5exon        model.(name).utr_exon] ;
  model.(link_name)(end+1,:) = [model.lengths.utr3exon        model.(name).utr_exon] ;
  if model.use.segments.polya_tail
    model.(link_name)(end+1,:) = [model.lengths.polya_tail      model.(name).utr_exon] ;
  end
  %model.(link_name)(end+1,:) = [model.lengths.intergenic_long model.(name).intergenic] ;
  %if model.use.segments.intercistronic
  %  model.(link_name)(end+1,:) = [model.lengths.intercistronic  model.(name).intercistronic] ;
  %end
  if isfield(model.use.segments, 'rna_seq_polya')&&model.use.segments.rna_seq_polya
    model.(link_name)(end+1,:) = [model.lengths.rna_seq_polya  model.(name).rna_seq_polya] ;
  end
if isfield(model.use, 'non_coding') && model.use.non_coding
	model.(link_name)(end+1,:) = [model.lengths.first_nc_exon  model.(name).utr_exon] ;
	model.(link_name)(end+1,:) = [model.lengths.middle_nc_exon  model.(name).utr_exon] ;
	model.(link_name)(end+1,:) = [model.lengths.last_nc_exon  model.(name).utr_exon] ;
end

end

%% mappings between length and intron list
for j=1:length(model.segment_feature_names)
  name = sprintf('segment_feature_%i',j);
  link_name = [name '_links'];
  model.(link_name) = [] ;   
  model.(link_name)(end+1,:) = [model.lengths.intron          model.(name).intron] ;
end
%% mappings between length and intron quality
for j=1:length(model.segment_feature_names)
  name = sprintf('segment_score_%i',j);
  link_name = [name '_links'];
  model.(link_name) = [] ;   
  model.(link_name)(end+1,:) = [model.lengths.intron          model.(name).intron] ;
end

%%mappings between length and segments

model.seg_links = [] ;
model.seg_links(end+1,:) = [model.lengths.intergenic      model.segments.intergenic] ;
model.seg_links(end+1,:) = [model.lengths.first_cds_exon  model.segments.cds_exon] ;
model.seg_links(end+1,:) = [model.lengths.last_cds_exon   model.segments.cds_exon] ;
model.seg_links(end+1,:) = [model.lengths.middle_cds_exon model.segments.cds_exon] ;
model.seg_links(end+1,:) = [model.lengths.single_cds_exon model.segments.cds_exon] ;
model.seg_links(end+1,:) = [model.lengths.intron          model.segments.intron] ;
if model.use.segments.transexon
  model.seg_links(end+1,:) = [model.lengths.transexon       model.segments.transexon] ;
end
model.seg_links(end+1,:) = [model.lengths.utr5exon        model.segments.utr5exon] ;
model.seg_links(end+1,:) = [model.lengths.utr3exon        model.segments.utr3exon] ;
if model.use.segments.polya_tail
  model.seg_links(end+1,:) = [model.lengths.polya_tail      model.segments.polya_tail] ;
end
if model.use.segments.intercistronic
  model.seg_links(end+1,:) = [model.lengths.intergenictrans  model.segments.intergenictrans] ;
  model.seg_links(end+1,:) = [model.lengths.intercistronic  model.segments.intercistronic] ;
end
%model.seg_links(end+1,:) = [model.lengths.intergenic_long model.segments.intergenic] ;

if isfield(model.use, 'non_coding') && model.use.non_coding
	model.seg_links(end+1,:) = [model.lengths.first_nc_exon  model.segments.nc_exon] ;
	model.seg_links(end+1,:) = [model.lengths.middle_nc_exon  model.segments.nc_exon] ;
	model.seg_links(end+1,:) = [model.lengths.last_nc_exon  model.segments.nc_exon] ;
end

if isfield(model.use.segments, 'rna_seq_polya')&&model.use.segments.rna_seq_polya
  model.seg_links(end+1,:) = [model.lengths.rna_seq_polya  model.segments.rna_seq_polya] ;
end

