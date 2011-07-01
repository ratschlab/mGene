function model = define_contents(model)
  
% declare content types

fid = 1;  
model.mod_words = [];
% the first contents consider all k-mers. the last three
% frame-content sensors only consider codons, with an offset from
% the end of the segment

content_start =  model.cnt_plifs+1;
cnt=0;

model.contents.intergenic = content_start+cnt ;cnt=cnt+1;
model.mod_words = [model.mod_words; 1 0] ;
model.contents.utr_exon   = content_start+cnt ;cnt=cnt+1;
model.mod_words = [model.mod_words; 1 0] ;
model.contents.cds_exon   = content_start+cnt ;cnt=cnt+1;
model.mod_words = [model.mod_words; 1 0] ;
model.contents.intron     = content_start+cnt ;cnt=cnt+1;
model.mod_words = [model.mod_words; 1 0] ;
if model.use.segments.intercistronic
  model.contents.intercistronic = content_start+cnt;cnt=cnt+1;
  model.mod_words = [model.mod_words; 1 0] ;
end
if model.use.contents.frame
  model.contents.cds_frame0   = content_start+cnt ;cnt=cnt+1;
  model.mod_words = [model.mod_words; 3 0] ;
  model.contents.cds_frame1   = content_start+cnt ;cnt=cnt+1;
  model.mod_words = [model.mod_words; 3 1] ;
  model.contents.cds_frame2   = content_start+cnt ;cnt=cnt+1;
  model.mod_words = [model.mod_words; 3 2] ;
end
if isfield(model.use.contents, 'rna_seq_polya')&&model.use.contents.rna_seq_polya
  model.contents.rna_seq_polya = content_start+cnt ;cnt=cnt+1;
  model.mod_words = [model.mod_words; 1 0] ;
end
if isfield(model.use, 'non_coding') && model.use.non_coding
	model.contents.nc_exon = content_start+cnt ;cnt=cnt+1;
	model.mod_words = [model.mod_words; 1 0] ;
end
model.cnt_contents = length(fieldnames(model.contents));

%hack
if model.cnt_contents<8
  fprintf(fid,'introducing dummy contents\n')
  for i= 1:8-model.cnt_contents
    if i==1
      model.contents.('dummy')= content_start+cnt;cnt=cnt+1;
    else
      model.contents.(['dummy' num2str(i)])= content_start+cnt;cnt=cnt+1;
    end
  end
  model.cnt_contents = length(fieldnames(model.contents));
  model.mod_words = [model.mod_words; 1 0] ;
end

% declare content range
model.contents_range = struct;
content_names = fieldnames(model.contents) ;
for s = 1:length(content_names)
  model.contents_range = setfield(model.contents_range, content_names{s},[-10,10]);
end

model.contents_fields = fieldnames(model.contents) ;

model.cnt_plifs = model.cnt_plifs + model.cnt_contents;
