function [xis,parameters,hvar,H,obj] = res2param(res, num_examples, model, QP, C, weight_names)
% [xis,parameters,hvar,H] = res2param(res, num_examples, model)
 

%%%%%% decompose res 

vector = res(1:model.cnt_parameters)' ;
parameters = init_weights(model) ;
param_name = weight_names ;%fieldnames(parameters);
last_pos = 0 ;
for p = 1:length(param_name)
  old_p = getfield(parameters, param_name{p}) ;
  parameters = setfield(parameters, param_name{p}, vector(last_pos+1:last_pos+length(old_p)));
  last_pos = last_pos + length(old_p) ;
end
assert(last_pos==length(vector)) ;


xis = res(model.cnt_parameters+1:model.cnt_parameters+num_examples);
 
hvar = res(model.cnt_parameters+num_examples+1:end);
H = [] ;
obj = struct;
%for i=1:length(model.plif2param{1})
%  H(i).hvar = hvar(model.plif2param{i});
%end


if ~isfield(QP, 'Q')
  return
end

%%%%%% decompose obj 

r = res;
r(1:model.cnt_parameters) = 0;
r(end-length(hvar)+1:end) = 0;
obj.slack = sum(r.*QP.f) ;

r = res;
r(model.cnt_parameters-model.cnt_transitions+1:end) = 0;
obj.plif_ys_sq = 0.5*r'*QP.Q*r;

r = res;
r(1:model.cnt_parameters-model.cnt_transitions) = 0;
r(model.cnt_parameters+1:end) = 0;
obj.transitions_sq = 0.5*r'*QP.Q*r;

r = res;
nof_smooth =(model.cnt_plifs*(model.bins-1))+model.cnt_transitions;% because of trans_smooth
r(1:end-nof_smooth) = 0;
obj.smoothness_sq = 0.5*r'*QP.Q*r;

r(end-model.cnt_transitions+1:end) = 0;
obj.smooth  = sum(r.*QP.f) ;

r = res;
r(1:model.cnt_parameters+num_examples+nof_smooth-model.cnt_transitions) = 0;
obj.transitions = QP.f'*r ;


%%%


obj.OBJ  = obj.slack + obj.smooth + obj.plif_ys_sq + obj.transitions_sq + obj.smoothness_sq + obj.transitions ;

OBJ_total = sum(res.*QP.f) + 0.5*res'*QP.Q*res ;

if (abs(OBJ_total-obj.OBJ)>1e-8)
  OBJ_total-obj.OBJ
  %keyboard
  obj.OBJ = OBJ_total 
end ;

%%%
