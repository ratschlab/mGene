function diff_fields = compare_PAR(PAR_1,PAR_2,silent)

if nargin<3
  silent=0;
end

diff_fields ={};
if isequal(PAR_1,PAR_2)
  fprintf(1,'\t parameters identical\n')
else
  fields_1 = fieldnames(PAR_1);
  fields_2 = fieldnames(PAR_2);
  common_fields = intersect(fields_1,fields_2);
  for i=1:length(common_fields)
    par_1 = getfield(PAR_1,common_fields{i});
    par_2 = getfield(PAR_2,common_fields{i});
    if ~isequal(par_1,par_2)
      diff_fields{end+1} = common_fields{i};
      if ~silent
        fprintf(1,'\t field "%s" not identical\n',common_fields{i})
        fprintf(1,'\t %s, 1:\n',common_fields{i});disp(par_1)
        fprintf(1,'\t %s, 2:\n',common_fields{i});disp(par_2)
      end
    end
  end
  fields_x1= setdiff(fields_1,fields_2);
  if ~isempty(fields_x1)
    if ~silent
      fprintf(1,'\t additional fields in PAR_1: \n')
    end
    diff_fields = {diff_fields{:},fields_x1{:}};
  end
  fields_x2 = setdiff(fields_2,fields_1);
  if ~isempty(fields_x2)
    if ~silent
      fprintf(1,'\t additional fields in PAR_2: \n')
    end
    diff_fields = {diff_fields{:},fields_x2{:}};
  end
end
   