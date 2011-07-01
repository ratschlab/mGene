function [model] = make_model(organism,LSL,Signals)
% [model] = make_model(bins)

model.bins = LSL.method.plif_bins ;
% model.use_gc =  LSL.method.use_gc;
model.use = LSL.method.use;
model.track_names = LSL.method.track_names;
model.track_monoton_functions = LSL.method.track_monoton_functions;
model.segment_feature_names = LSL.method.segment_feature_names;
model.segment_feature_monoton_functions = LSL.method.segment_feature_monoton_functions;

polya_sixmers = Signals.polya.consensus;

model = define_segments(model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%define features

model.cnt_plifs = 0;
model = define_contents(model);
model = define_lengths(model,organism);
model = define_signals(model,polya_sixmers);
model = define_addContentFeatures(model);
model = define_intron_lists(model);
model = define_monotonicity_constraints(model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%define states

model = define_states(model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define mappings between length, contents and segments

model = define_mappings(model);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = gen_transition_matrix(model) ;
% transition_pointers used to be called penalties
model = define_transition_pointers(model);
model.A(~isinf(model.A))=0 ;



model.cnt_parameters = model.cnt_transitions + model.cnt_plifs*model.bins ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define loss 
model =  define_loss(model);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mappings of PLIFS to parameters

model = plif_to_param_mapping(model);
