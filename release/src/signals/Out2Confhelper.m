function [dummy1, dummy2]=Out2Confhelper(PAR)
% Out2Confhelper(PAR)
try
  pause('on')
  warning('off', 'MATLAB:typeaheadBufferOverflow');  
catch
  % works only for matlab
end

dummy1=[] ;
dummy2=[] ;

try
  load(PAR.fn_pred,'output') ;
catch
  pause(60)
  load(PAR.fn_pred,'output') ;
end
disp('Convert SVM Output to Confidence on fn_pred')
%[Conf, Conf_cum] = Out2Conf(output, PAR.prob, PAR.prob_cum, PAR.limits);

if ~isempty(output)
  limits = PAR.limits(~isinf(PAR.limits));
  prob = PAR.prob(~isinf(PAR.limits));
  Conf = plif_transform(output, horz(prob), horz(limits));
  
  prob_cum = PAR.prob_cum(~isinf(PAR.limits));
  Conf_cum = plif_transform(output, horz(prob_cum), horz(limits));
else
  Conf = [];
  Conf_cum = [];
end
if ~(all(Conf>=0-1e-8&Conf<=1+1e-8))
  save('~/tmp/ws_Out2Confhelper_error');
  error('Conf out of range [0, 1]');
end
if ~(all(Conf_cum>=0-1e-8&Conf_cum<=1+1e-8))
  save('~/tmp/ws_Out2Confhelper_error');
  error('Conf_cum out of range [0, 1]');
end

fprintf('append Conf and Conf_cum to file %s ', PAR.fn_pred)
%save(PAR.fn_pred,'-append','Conf','Conf_cum');
inventory = {'output', 'pos', 'Conf', 'Conf_cum'};
save_append(PAR.fn_pred, 1, 'inventory', inventory,'Conf', Conf, 'Conf_cum', Conf_cum);

return

function vec = horz(vec)

if size(vec, 1)==1
  return;
elseif size(vec, 2)==1
  vec = vec';
else
  error('row or column vector expected: but size is: [%ix%i]', size(vec, 1), size(vec, 2));
end

return
