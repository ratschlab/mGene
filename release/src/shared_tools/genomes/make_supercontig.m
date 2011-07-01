function [] = make_supercontig(pattern, list, save_name)
%function [] = make_supercontig(pattern, list, save_name)
% <pattern> has to be a string 'pre*post' where '*' is the wild card sign
%
% <list> is a list of strings ordered in the way the they should exchange 
%  the wild card in pattern.
%
% <save_name> is the file where the supercontig has to be written to
  if nargin < 3
    error('not enough parameters')
  end
  wc=find(pattern=='*');
  if isempty(wc)
    error('No wild-card character found.');
  end
  pre=pattern(1:wc-1);
  post=pattern(wc+1:end);
  cmd_str='cat ';
  for i=1:length(list)
    ext_str=sprintf('%s%s%s ',pre,list{i},post);
    cmd_str=[cmd_str,ext_str];
  end
  cmd_str=[cmd_str,' > ',save_name];
  unix(cmd_str);
end