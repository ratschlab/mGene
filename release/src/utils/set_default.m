function new_field = set_default(old_struct,field,defaults,silent)  ;

% new_struct = set_default(old_struct,defaults,overwrite)  ;
% sets all values of all fields in old_struct to the default values specified in
% defaults IF field is empty in old_struct. If overwrite is set all entries
% will be overwritten.
% defaults must be a struct.
%
% e.g. 
% default.use_cons = 0;
% default.plif.bins = 50;
% method.use_cons = 0;
% method = set_default(method,default,1) 

if ~isfield(old_struct,field)
  new_field = defaults;
  return
end  

old_struct = old_struct.(field);

if nargin<4
  silent = 1 ;
end

names = fieldnames(defaults);
for j=1:length(names)
  if ~isfield(old_struct,names{j}) || isempty(old_struct.(names{j}))||...
        (iscell(old_struct.(names{j}))&&length(old_struct.(names{j}))==1&&isempty(old_struct.(names{j}){1}))
    if ~isfield(old_struct,names{j}) && ~silent
      warning('new field added')
      names{j}
    end
    
    old_struct.(names{j}) = defaults.(names{j});
  elseif isstruct(defaults.(names{j}))
    old_struct.(names{j}) = set_default(old_struct,names{j},defaults.(names{j}));
    %names2 =  fieldnames(defaults.(names{j}));
    %for k= 1:length(names2)
    %  if ~isfield(old_struct.(names{j}),names2{k})||isempty(old_struct.(names{j}).(names2{k}))||overwrite
    %    old_struct.(names{j}).(names2{k}) = defaults.(names{j}).(names2{k});
    %  end
    %end
  end
end
% if length(fieldnames(PAR))>length(names)
%   keyboard
% end
% PAR = orderfields(PAR,defaults);

new_field = orderfields(old_struct,defaults);
