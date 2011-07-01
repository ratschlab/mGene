function model=fix_model(model) ;
% model = fix_model(model) ;

model.contents = reorder_fields(model.contents, model.contents_fields) ;
model.lengths = reorder_fields(model.lengths, model.lengths_fields) ;
model.signals = reorder_fields(model.signals, model.signals_fields) ;
model.segments = reorder_fields(model.segments, model.segments_fields) ;
