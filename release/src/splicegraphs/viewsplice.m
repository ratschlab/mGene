function viewsplice(truepos, name, labels,colors, strands)

start=inf ;
stop=-inf ;
for k=1:length(truepos)
  start=min(start,truepos{k}(1,1)-50) ;
  stop=max(stop,truepos{k}(end,2)+50) ;
end ;

cnt=0;
width=0.075 ; width2=0.03; step=0.75 ;
for i=ceil(start/1000)*1000:1000:floor(stop/1000)*1000,
  plot([i i], [step/4 step*length(truepos)+step/2], 'k:') ;
  hold on
  h=text(i, step/8, sprintf('%i', i)) ;
  set(h, 'HorizontalAlignment', 'center') ;
  set(h, 'FontSize', 7) ;
end ;

if nargin<4 | isempty(colors)
  colors={} ;
  for i=1:500,
    colors{i,1}='k' ;
    colors{i,2}=[0 0 0] ;
  end ;
end ;
if nargin<3,
  labels={'Annotation', 'EST', 'Prediction'} ;
end ;

h=text(start, cnt+(length(truepos)+1)*step, name) ;
set(h, 'HorizontalAlignment', 'center') ;
set(h, 'FontSize', 10) ;
for k=1:length(truepos)
  
  cnt=cnt+step ;
  
  plot([start truepos{k}(1,1)], [cnt cnt], colors{k,1}) ;
  hold on
  plot([truepos{k}(end,2) stop], [cnt cnt], colors{k,1}) ;
  axis off
  h=text(start, cnt, [labels{k} ' ']) ;
%  set(h, 'VerticalAlignment', 'center') ;
  set(h, 'HorizontalAlignment', 'right') ;
  set(h, 'FontSize', 7) ;
  h=text(stop, cnt, ['  ' strands(k)]) ;
%  set(h, 'VerticalAlignment', 'center') ;
  set(h, 'HorizontalAlignment', 'left') ;
  set(h, 'FontSize', 7) ;
  axis([start stop 0 cnt+step]) ;
  
  for i=1:size(truepos{k},1)
    if size(truepos{k},2)==2 | truepos{k}(i,3)==1
      patch([truepos{k}(i,1) truepos{k}(i,1) truepos{k}(i,2) truepos{k}(i,2)],...
            [cnt-width, cnt+width, cnt+width, cnt-width], colors{k,1}) ;
    else
      patch([truepos{k}(i,1) truepos{k}(i,1) truepos{k}(i,2) truepos{k}(i,2)],...
            [cnt-width2, cnt+width2, cnt+width2, cnt-width2], colors{k,1}) ;
    end ;
    h=text(mean([truepos{k}(i,1) truepos{k}(i,2)]), cnt-3*width, sprintf('%i', truepos{k}(i,2)-truepos{k}(i,1))) ;
    set(h, 'HorizontalAlignment', 'center') ;
    set(h, 'FontSize', 7) ;
  end ;
  
  for i=1:size(truepos{k},1)-1
    h=line([truepos{k}(i,2) mean([truepos{k}(i,2) truepos{k}(i+1,1)]) truepos{k}(i+1,1)],...
         [cnt+width, cnt+2*width, cnt+width]) ;
    set(h,'color',colors{k,2}) ;
    h=text(mean([truepos{k}(i,2) truepos{k}(i+1,1)]), cnt+4*width, sprintf('%i', truepos{k}(i+1,1)-truepos{k}(i,2))) ;
    set(h, 'HorizontalAlignment', 'center') ;
    set(h, 'FontSize', 7) ;
  end ;
end ;
  
return

close all
names={'F49A5.7', 'C05A9.3', 'F26D2.12', 'F40G9.5', ...
           'M03E7.3', 'T13B5.7', 'F40H7.5', 'K08D10.9', 'T13B5.7', 'M4.1', ...
           'Y69H2.7', 'C35E7.7', 'F21C10.9', 'B0432.7', ...
           'T04C10.3', 'F07E5.6', 'K04E7.1', 'C18F10.9', 'F59H6.1', 'C18G1.8',...
           'F49H6.10', 'T12C9.7', 'F02D10.3'} ;
q=0; 
for i=1:4:23,
  figure
  q=q+1 ;
  for j=i:i+3,
    if j<=23,
      subplot(4,1,j-i+1) ;
      viewsplice({anno{j}, est{j}, pred{j}},names{j}) ;
      
    end ;
  end ;
  pause(3) ;
  print('-depsc2', sprintf('splice%i.eps', q)) ;
end ;
