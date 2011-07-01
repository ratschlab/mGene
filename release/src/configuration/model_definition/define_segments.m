function model = define_segments(model)
% declare segment types

  cnt=0;  
model.segments.intergenic = cnt ;cnt=cnt+1;
model.segments.utr5exon   = cnt ;cnt=cnt+1;
if model.use.segments.transexon 
  model.segments.transexon  = cnt ;cnt=cnt+1;
end
model.segments.cds_exon   = cnt ;cnt=cnt+1;
model.segments.utr3exon   = cnt ;cnt=cnt+1;
if model.use.segments.polya_tail
  model.segments.polya_tail = cnt ;cnt=cnt+1;
end
if model.use.segments.intercistronic
  model.segments.intercistronic  = cnt ;cnt=cnt+1;
  model.segments.intergenictrans  = cnt ;cnt=cnt+1;
end
model.segments.intron     = cnt ;cnt=cnt+1;
if isfield(model.use.segments, 'rna_seq_polya')&&model.use.segments.rna_seq_polya
  model.segments.rna_seq_polya = cnt ;cnt=cnt+1;
end

if isfield(model.use, 'non_coding') && model.use.non_coding
	model.segments.nc_exon = cnt ;cnt=cnt+1;
end

model.cnt_segments = length(fieldnames(model.segments));
model.segments_fields = fieldnames(model.segments) ;
