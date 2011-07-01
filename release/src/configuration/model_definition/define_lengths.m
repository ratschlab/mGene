function model = define_lengths(model,organism)
  
fid = 1; 
% declare length types
length_start = model.cnt_plifs+1 ;
cnt = 0;

model.lengths.intergenic = length_start+cnt; cnt = cnt+1;
model.lengths.utr5exon = length_start+cnt ;cnt = cnt+1;
if model.use.segments.transexon 
  model.lengths.transexon = length_start+cnt ;cnt = cnt+1;
end
if isfield(model.use, 'non_coding') && model.use.non_coding
	model.lengths.first_nc_exon = length_start+cnt ;cnt = cnt+1;
	model.lengths.middle_nc_exon = length_start+cnt ;cnt = cnt+1;
	model.lengths.last_nc_exon = length_start+cnt ;cnt = cnt+1;
end
model.lengths.single_cds_exon = length_start+cnt ;cnt = cnt+1;
model.lengths.first_cds_exon = length_start+cnt ;cnt = cnt+1;
model.lengths.middle_cds_exon = length_start+cnt ;cnt = cnt+1;
model.lengths.last_cds_exon = length_start+cnt ;cnt = cnt+1;
model.lengths.utr3exon = length_start+cnt;cnt = cnt+1;
if model.use.segments.polya_tail 
  model.lengths.polya_tail = length_start+cnt;cnt = cnt+1;
end
model.lengths.intron = length_start+cnt ;cnt = cnt+1;
if model.use.segments.intercistronic
  model.lengths.intercistronic = length_start+cnt ;cnt = cnt+1;
  model.lengths.intergenictrans = length_start+cnt ;cnt = cnt+1;
end
model.lengths.intergenic_long = length_start+cnt ;cnt = cnt+1;
if isfield(model.use.segments, 'rna_seq_polya')&&model.use.segments.rna_seq_polya
  model.lengths.rna_seq_polya = length_start+cnt ;cnt = cnt+1; 
end

model.cnt_lengths = length(fieldnames(model.lengths));


% declare length range
s1=which([organism.name '_lengths']);
s2=which([organism.clade '_lengths']);
if ~isempty(s1)
  fprintf(fid,'taking lengths range from file: %s_lengths.m \n',organism.name)
  model.lengths_range = feval([organism.name '_lengths']);
elseif ~isempty(s2)
  fprintf(fid,'taking lengths range from file: %s_lengths.m \n',organism.clade)
  model.lengths_range = feval([organism.clade '_lengths']);
else
  fprintf(fid,'no lengths range specified, using default settings\n')
  model.lengths_range = feval('default_lengths');
end
model.lengths_fields = fieldnames(model.lengths) ;

model.cnt_plifs = model.cnt_plifs + model.cnt_lengths;

