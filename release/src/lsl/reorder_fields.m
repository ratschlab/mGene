function n = reorder_fields(s, fields)
% n = reorder_fields(s, fields)

  if isempty(s),
    n=s ;
    return ;
  end ;

  n=struct ;
  for i=1:length(fields)
    for j=1:length(s) 
      n(j).(fields{i}) = s(j).(fields{i}) ;
    end ;
  end ;

  remaining = setdiff(fieldnames(s), fields) ;

  if ~isempty(remaining)
    remaining
    error('these fields are unassigned');
  end ;

