function [attributes, type_name, elems, action] = parse_starter(fd, types, types_IDs, contig, source)
% parse_starter(fd)

attributes = {} ;
type_name = '' ;
elems = {} ;
action = '' ;


Line = fgetl(fd);
elems = separate(Line);
assert(~isempty(elems{1}));

if Line(1)=='#',
  action = 'continue' ;
  return ;
end
if ~isempty(contig)
  if ~isequal(upper(elems{1}),upper(contig)),
    action='continue' ;
    return ;
  end
end
if ~ischar(Line),
  action='break' ;
  return ;
end

if isempty(elems{4}) | isempty(elems{4})
  action = 'continue' ;
  return ;
end


assert(str2double(elems{4})<=str2double(elems{5}) || str2double(elems{5})<0)

if~(elems{7}=='+' || elems{7}=='-') ;
  if elems{7}~='.',
    warning('strand missing')
  end ;
end

if ~isempty(source)
  if ~isequal(elems{2},source)
    error('something wrong with source')
  end  
end

if ~ismember(elems{3},types) && ~ismember(elems{3},types_IDs)
  action = 'continue' ;
  return ;
end
attributes = separate(elems{9}, ';');

if ismember(elems{3},types);
  type_name = elems{3};
elseif ismember(elems{3},types_IDs);
  type_name =types(isequal(elems{3},types_IDs));
else
  error('wrong type')
end

action='' ;
