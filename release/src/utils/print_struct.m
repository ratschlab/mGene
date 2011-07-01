function print_struct(fid,s)
% print_struct(fid,s)

fields = fieldnames(s);
max_len = 0;
for f=1:length(fields)
  max_len = max(max_len, length(fields{f}));
end

for f=1:length(fields),
  value = getfield(s,fields{f});
  if isempty(value),
    fprintf(fid, '%*s: %s\n',max_len+3,fields{f},'[]');
  elseif iscell(value),
    if size(value,1)==1 & size(value,2)<10,
      if isscalar(value{1}),
        vstr = ['{[' num2str(value{1}) ']'];
        for i=2:size(value,2),
          vstr = [vstr ' [' num2str(value{i}) ']'];
        end
        vstr = [vstr '}'];
      else
        if isempty(value{1}),
          vstr = ['{[]'];
        else
          vstr = ['{[' class(value{1}) ']'];
        end
        for i=2:size(value,2),
          if isempty(value{i}),
            vstr = [vstr ' []'];
          else
            vstr = [vstr ' [' class(value{i}) ']'];
          end
        end
        vstr = [vstr '}'];        
      end
      fprintf(fid, '%*s: %s\n',max_len+3,fields{f},vstr);  
    else
      fprintf(fid, '%*s: %s\n',max_len+3,fields{f},...
              sprintf('{%ix%i cell}',size(value,1),size(value,2)));  
    end
  elseif isstruct(value),
    fprintf(fid, '%*s: %s\n',max_len+3,fields{f},...
            sprintf('[%ix%i struct]',size(value,1),size(value,2)));
  elseif isscalar(value),
    if ischar(value),
      fprintf(fid, '%*s: "%s"\n',max_len+3,fields{f}, value);
    elseif isinteger(value) | mod(value,1)==0,
      fprintf(fid, '%*s: %i\n',max_len+3,fields{f}, value);
    elseif isnumeric(value),
      fprintf(fid, '%*s: %f\n',max_len+3,fields{f}, value);
    else 
      fprintf(fid, '%*s: ...\n',max_len+3,fields{f});     
    end    
  elseif isvector(value),
    if size(value,1)==1 & size(value,2)<10,
      if ischar(value),
        fprintf(fid, '%*s: "%s"\n',max_len+3,fields{f}, value);
      elseif isscalar(value(1)),
        vstr = ['[' num2str(value(1))];
        for i=2:size(value,2),
          vstr = [vstr ' ' num2str(value(i))];
        end
        vstr = [vstr ']'];
      else
        vstr = ['[' class(value(1))];
        for i=2:size(value,2),
          vstr = [vstr ' ' class(value(i))];
        end
        vstr = [vstr ']'];        
      end
      fprintf(fid, '%*s: %s\n',max_len+3,fields{f},vstr);      
    else
      fprintf(fid, '%*s: %s\n',max_len+3,fields{f},...
              sprintf('[%ix%i %s]',size(value,1),size(value,2),class(value)));
    end
  else
    fprintf(fid, '%*s: %s\n',max_len+3,fields{f},...
              sprintf('[%ix%i %s]',size(value,1),size(value,2),class(value)));      
  end
end