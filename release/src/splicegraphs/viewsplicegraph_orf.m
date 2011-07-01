function viewsplicegraph_orf(gene,flag,param)

% viewsplicegraph_orf(gene,flag,param)
% set flag = 'g' to view the graph
% set flag = 't' to view the transcripts
% set flag = 'o' to view the open reading frames
%
% Appending flags views all, e.g. 'go' views both graphs and open reading frames.
% Default behaviour (no flag) views everything.
%
% param is an optional structure with fields and default values:
%   FontSize     8
%   FontWeight   'normal'
%   LineWidth    1
%   grid_on      1
%   colors       {[0 0 0]   transcripts exon 
%                 [0 0 0]   transcripts intron
%                 [1 0 0]   summary exon
%                 [0 0 0]}  summary intron
%
% modified from cheng's version by GS 17/08/06
% => intron retention represented with both exons on the same level    
% => param added 
%
% modified by Cheng Soon Ong on 17/04/07
% - split the viewing code into sub functions.
% - display also the open reading frames

if nargin<3,param=struct; end ;
if nargin<2,flag='gto' ; end ;



%%%%% SET PARAMETERS
%if length(strands)<length(labels)
%  strands = repmat(strands(1),1,length(labels)) ;
%end
  

if ~exist('param')
  param = struct;
end
param.width=0.075;
param.width2=0.03;
param.step=0.4;
param.padlength = 250;

if ~isfield(param,'flip')
  param.flip = 0;
end
if ~isfield(param,'FontSize')
  param.FontSize = 8;
end
if ~isfield(param,'FontWeight')
  param.FontWeight = 'normal';
end
if ~isfield(param,'LineWidth')
  param.LineWidth = 1;
end 
if ~isfield(param,'grid_on')
  param.grid_on = 1;
end
if ~isfield(param,'colors')
  param.colors={} ;
  for i=1:500,
    param.colors{i,1}=[0 0.3 0] ;% transcripts exon 
    param.colors{i,2}=[0 0.3 0] ;% transcripts intron
    param.colors{i,3}=[1 0 0] ;% summary exon
    param.colors{i,4}=[0 0 0] ;% summary intron
    param.colors{i,5}=[0 0 1] ;% orf exon
    param.colors{i,6}=[0 0 1] ;% orf intron
  end 
end
if ~isfield(param,'maxdisplay')
  param.maxdisplay = 200;
end


%%%%% DETERMINE BOUNDARIES
start=inf ;
stop=-inf ;
for k=1:length(gene.exons)
  start=min(start,gene.exons{k}(1,1)-param.padlength) ;
  stop=max(stop,gene.exons{k}(end,2)+param.padlength) ;
end ;

for k=1:size(gene.splicegraph{1},1)
  start=min(start,min(gene.splicegraph{1}(1,:)-param.padlength)) ;
  stop=max(stop,max(gene.splicegraph{1}(2,:)+param.padlength)) ;
end ;
param.start = start;
param.stop = stop;

for i=1:length(gene.transcripts)
  gene.transcripts{i} = strrep(gene.transcripts{i}, '_', '\_') ;
  idx=find(gene.transcripts{i}=='-') ;
  if length(idx)>=2,
    gene.transcripts{i} = gene.transcripts{i}(idx(1)+1:idx(end)-1) ;
  end ;
  idx=find(gene.transcripts{i}=='|') ;
  if length(idx)>=2,
    gene.transcripts{i} = gene.transcripts{i}(idx(1)+1:idx(end)-1) ;
  end ;
end 

clf;hold on
axis off
cnt = 0;
if strfind(flag,'t')
  cnt = viewtranscripts(gene.exons,gene.exons_confirmed,gene.transcripts,gene.strands,param,cnt);
end
if strfind(flag,'o')
  cnt = vieworf(gene.splicegraph{1},gene.splicegraph{2},gene.chr,gene.strands(1),param,cnt);
end
if strfind(flag,'g')
  cnt = viewgraph(gene.splicegraph{1},gene.splicegraph{2},gene.chr,gene.strands(1),param,cnt);
end
cnt = viewname(gene.name,gene.chr,gene.strands(1),param,cnt);



return




function cnt = viewtranscripts(truepos,confirmed,labels,strands,param,cnt)
% View the individual transcripts that were used to infer the splice graphs

% exons_confirmed{idx1}(idx2,1) gives the number of exon confirmations
% exons_confirmed{idx1}(idx2,2) gives the number of intron confirmations

start = param.start;
stop = param.stop;
width = param.width;
width2 = param.width2;
step = param.step;

if strands(1) == '-' && param.flip
  for k=1:length(truepos)
    truepos{k} = -truepos{k};
  end
end

if length(truepos)>param.maxdisplay,
  cnt = cnt + step ;
  h = text(mean([start stop]), cnt+0*width, sprintf('%i more transcripts not displayed', length(truepos)-param.maxdisplay), 'HorizontalAlignment', 'center', 'FontSize',param.FontSize, 'FontWeight',param.FontWeight) ;
end ;

rand('seed', 34834) ;
pp=randperm(length(truepos)) ;

for k=pp(1:min(param.maxdisplay,length(truepos)))
  cnt=cnt+step ;
    
  %%%%%% DRAW LINE FROM GENE START TO START OF FIRST EXON 
  if strands(1) == '+' || ~param.flip
    h = plot([start truepos{k}(1,1)], [cnt cnt], 'k', 'LineWidth',0.5) ;
  else
    h = plot([-start truepos{k}(1,1)], [cnt cnt], 'k', 'LineWidth',0.5) ;
  end    
  %set(h, 'LineWidth',0.5) ;
  hold on
  %%%%%% DRAW LINE FROM END OF LAST EXON TO GENE END 
  if strands(1) == '+' || ~param.flip
    h = plot([truepos{k}(end,2) stop], [cnt cnt],'k', 'LineWidth', 0.5) ;
  else
    h = plot([truepos{k}(end,2) -stop], [cnt cnt],'k', 'LineWidth', 0.5) ;
  end
  %  set(h, 'LineWidth',0.5) ;
  axis off
  %%%%%% ADD TRANSCRIPT NAME 
  if strands(1) == '+' || ~param.flip
    h = text(start, cnt, [labels{k} ' ']) ;
  else
    h = text(-stop, cnt, [labels{k} ' ']) ;
  end
  set(h, 'HorizontalAlignment', 'right') ;
  set(h, 'FontSize',param.FontSize ) ;
  set(h, 'FontWeight',param.FontWeight) ;
  %%%%%% ADD STRAND 
  if strands(1) == '+' || ~param.flip
    h = text(stop, cnt, ['  ' strands(k)]) ;
  else
    h = text(-start, cnt, ['  ' strands(k)]) ;
  end
  set(h, 'HorizontalAlignment', 'left') ;
  set(h, 'FontSize',param.FontSize ) ;
  set(h, 'FontWeight',param.FontWeight) ;
  
  %%%%%% DRAW EXONS 
  for i=1:size(truepos{k},1)
    if size(truepos{k},2)==2 | truepos{k}(i,3)==1
      patch([truepos{k}(i,1) truepos{k}(i,1) truepos{k}(i,2) truepos{k}(i,2)],...
	    [cnt-width, cnt+width, cnt+width, cnt-width],param.colors{mod(k,499)+1,1}) ;
    else
      patch([truepos{k}(i,1) truepos{k}(i,1) truepos{k}(i,2) truepos{k}(i,2)],...
	    [cnt-width2, cnt+width2, cnt+width2, cnt-width2],param.colors{mod(k,499)+1,1}) ;
    end ;
    %%%%%% ADD NUMBER OF CONFIRMATIONS
    h = text(mean([truepos{k}(i,1) truepos{k}(i,2)]), cnt-2*width, sprintf('%i', confirmed{k}(i,1))) ;
    set(h, 'HorizontalAlignment', 'center') ;
    set(h, 'FontSize',param.FontSize) ;
    set(h, 'FontWeight',param.FontWeight) ;
  end ;
  
  %%%%%% DRAW INTRONS 
  for i=1:size(truepos{k},1)-1
    h = line([truepos{k}(i,2) mean([truepos{k}(i,2) truepos{k}(i+1,1)]) truepos{k}(i+1,1)],...
             [cnt+width/2, cnt+2*width, cnt+width/2], 'color',param.colors{mod(k,499)+1,2}, 'LineWidth',param.LineWidth) ;
    %set(h,'color',param.colors{mod(k,499)+1,2}) ;
    %set(h, 'LineWidth',param.LineWidth) ;
    %%%%%% ADD NUMBER OF CONFIRMATIONS
    h = text(mean([truepos{k}(i,2) truepos{k}(i+1,1)]), cnt+0*width, sprintf('%i', confirmed{k}(i,2)), 'HorizontalAlignment', 'center', 'FontSize',param.FontSize, 'FontWeight',param.FontWeight) ;
    %set(h, 'HorizontalAlignment', 'center') ;
    %set(h, 'FontSize',param.FontSize) ;
    %set(h, 'FontWeight',param.FontWeight) ;
  end ;
end ;
%axis([start stop 0 cnt+step]) ;
  
function  cnt = vieworf(vertices,edges,chr,strand,param,cnt)
% View the top 5 open reading frames

param.num_orf = 5;
start = param.start;
stop = param.stop;
width = param.width;
width2 = param.width2;
step = param.step;

isoforms = find_all_isoforms(vertices,edges,strand);
[origisoforms,isoforms,isoORFlength,isostartpos,isostoppos] = find_all_orfs(isoforms,chr,strand,vertices);

cnt=cnt+step ;

for ix = min(param.num_orf,length(isoforms)):-1:1
  cnt=cnt+step ;
  to_draw = 0;
  %fprintf('cnt:%f\tisoform:%d\tstart:%d\tstop:%d\npath:',cnt,ix,isostartpos(ix),isostoppos(ix));
  %fprintf('%d,',isoforms{ix});
  %fprintf('\n');

  % draw the whole isoform
  for ixe = 1:length(origisoforms{ix})
    cur_vertex = vertices(:,origisoforms{ix}(ixe)) - [start;start];
    %%%%%% DRAW EXONS
    if strand == '+' || ~param.flip
      patch([cur_vertex(1)+start cur_vertex(1)+start,...
	     cur_vertex(2)+start cur_vertex(2)+start],...
	    [cnt-width, cnt+width, cnt+width, cnt-width],param.colors{1,1} );
    else
      patch(-[cur_vertex(1)+start cur_vertex(1)+start,...
	     cur_vertex(2)+start cur_vertex(2)+start],...
	    [cnt-width, cnt+width, cnt+width, cnt-width],param.colors{1,1} );
    end

    if ixe > 1
      %%%%%% DRAW INTRONS 
      istart = vertices(2,origisoforms{ix}(ixe-1));
      istop =  vertices(1,origisoforms{ix}(ixe));
      if strand == '+' || ~param.flip
	h = line([istart, mean([istart,istop]),istop],...
		 [cnt+width2,cnt+step*0.75,cnt+width2]);
      else
	h = line(-[istart, mean([istart,istop]),istop],...
		 [cnt+width2,cnt+step*0.75,cnt+width2], 'LineWidth',param.LineWidth, 'Color',param.colors{1,1});
      end
      %set(h, 'LineWidth',param.LineWidth) ;
      %set(h, 'Color',param.colors{1,1}) ;
    end
  end
  
  %%%%%% PRINT ORF LENGTH
  if strand == '+' || ~param.flip
    h = text(start-param.padlength, cnt, sprintf('%i ', isoORFlength(ix)));
  else
    h = text(-stop-param.padlength, cnt, sprintf('%i ', isoORFlength(ix)));
  end
  set(h, 'HorizontalAlignment', 'center') ;
  set(h, 'FontSize',param.FontSize) ;
  set(h, 'FontWeight',param.FontWeight) ;

  % draw only the open reading frame
  if length(isoforms{ix}) == 1
    %%%%%% DRAW EXONS
    if strand == '+' || ~param.flip
      patch([isostartpos(ix) isostartpos(ix),...
	     isostoppos(ix) isostoppos(ix)],...
	    [cnt-width, cnt+width, cnt+width, cnt-width],param.colors{1,5} );
    else
      patch(-[isostartpos(ix) isostartpos(ix),...
	     isostoppos(ix) isostoppos(ix)],...
	    [cnt-width, cnt+width, cnt+width, cnt-width],param.colors{1,5} );
    end
    %%%%%% ADD EXON LENGTH
    if strand == '+' || ~param.flip
      h = text(mean([isostoppos(ix) isostartpos(ix)]), cnt-2*width, ...
	       sprintf('%i', isostoppos(ix)-isostartpos(ix)+1));
    else
      h = text(-mean([isostoppos(ix) isostartpos(ix)]), cnt-2*width, ...
	       sprintf('%i', isostoppos(ix)-isostartpos(ix)+1));
    end
    set(h, 'HorizontalAlignment', 'center') ;
    set(h, 'FontSize',param.FontSize) ;
    set(h, 'FontWeight',param.FontWeight) ;
  else
    for ixe = 1:length(isoforms{ix})
      cur_vertex = vertices(:,isoforms{ix}(ixe)) - [start;start];
      if (cur_vertex(1)+start <= isostoppos(ix)) && (cur_vertex(2)+start >= isostoppos(ix))
	to_draw = 0;
	%%%%%% DRAW EXONS
	if strand == '+' || ~param.flip
	  patch([cur_vertex(1)+start cur_vertex(1)+start,...
		 isostoppos(ix) isostoppos(ix)],...
		[cnt-width, cnt+width, cnt+width, cnt-width],param.colors{1,5} );
	else
	  patch(-[cur_vertex(1)+start cur_vertex(1)+start,...
		 isostoppos(ix) isostoppos(ix)],...
		[cnt-width, cnt+width, cnt+width, cnt-width],param.colors{1,5} );
	end
	%%%%%% ADD EXON LENGTH
	if strand == '+' || ~param.flip
	  h = text(mean([isostoppos(ix) cur_vertex(1)+start]), cnt-2*width, ...
		   sprintf('%i', isostoppos(ix)-cur_vertex(1)-start+1));
	else
	  h = text(-mean([isostoppos(ix) cur_vertex(1)+start]), cnt-2*width, ...
		   sprintf('%i', isostoppos(ix)-cur_vertex(1)-start+1));
	end
	set(h, 'HorizontalAlignment', 'center') ;
	set(h, 'FontSize',param.FontSize) ;
	set(h, 'FontWeight',param.FontWeight) ;

	%%%%%% DRAW INTRONS 
	istart = vertices(2,isoforms{ix}(ixe-1));
	istop =  vertices(1,isoforms{ix}(ixe));
	if strand == '+' || ~param.flip
	  h = line([istart, mean([istart,istop]),istop],...
		   [cnt+width2,cnt+step*0.75,cnt+width2]);
	else
	  h = line(-[istart, mean([istart,istop]),istop],...
		   [cnt+width2,cnt+step*0.75,cnt+width2]);
	end	  
	set(h, 'LineWidth',param.LineWidth) ;
	set(h, 'Color',param.colors{1,6}) ;
	%%%%%% ADD INTRON LENGTH
	if strand == '+' || ~param.flip
	  h = text(mean([istart istop]), cnt+0*width, sprintf('%i', istop-istart));
	else
	  h = text(-mean([istart istop]), cnt+0*width, sprintf('%i', istop-istart));
	end
	set(h, 'HorizontalAlignment', 'center') ;
	set(h, 'FontSize',param.FontSize) ;
	set(h, 'FontWeight',param.FontWeight) ;
      end
      if to_draw
	%%%%%% DRAW EXONS
	if strand == '+' || ~param.flip
	  patch([cur_vertex(1)+start cur_vertex(1)+start,...
		 cur_vertex(2)+start cur_vertex(2)+start],...
		[cnt-width, cnt+width, cnt+width, cnt-width],param.colors{1,5} );
	else
	  patch(-[cur_vertex(1)+start cur_vertex(1)+start,...
		 cur_vertex(2)+start cur_vertex(2)+start],...
		[cnt-width, cnt+width, cnt+width, cnt-width],param.colors{1,5} );
	end
	%%%%%% ADD EXON LENGTH
	if strand == '+' || ~param.flip
	  h = text(mean([cur_vertex(1)+start cur_vertex(2)+start]), cnt-2*width, ...
		   sprintf('%i', cur_vertex(2)-cur_vertex(1)+1));
	else
	  h = text(-mean([cur_vertex(1)+start cur_vertex(2)+start]), cnt-2*width, ...
		   sprintf('%i', cur_vertex(2)-cur_vertex(1)+1));
	end
	set(h, 'HorizontalAlignment', 'center') ;
	set(h, 'FontSize',param.FontSize) ;
	set(h, 'FontWeight',param.FontWeight) ;

	%%%%%% DRAW INTRONS 
	istart = vertices(2,isoforms{ix}(ixe-1));
	istop =  vertices(1,isoforms{ix}(ixe));
	if strand == '+' || ~param.flip
	  h = line([istart, mean([istart,istop]),istop],...
		   [cnt+width2,cnt+step*0.75,cnt+width2]);
	else
	  h = line(-[istart, mean([istart,istop]),istop],...
		   [cnt+width2,cnt+step*0.75,cnt+width2]);
	end
	set(h, 'LineWidth',param.LineWidth) ;
	set(h, 'Color',param.colors{1,6}) ;
	%%%%%% ADD INTRON LENGTH
	if strand == '+' || ~param.flip
	  h = text(mean([istart istop]), cnt+0*width, sprintf('%i', istop-istart));
	else
	  h = text(-mean([istart istop]), cnt+0*width, sprintf('%i', istop-istart));
	end
	set(h, 'HorizontalAlignment', 'center') ;
	set(h, 'FontSize',param.FontSize) ;
	set(h, 'FontWeight',param.FontWeight) ;
      end
      if (cur_vertex(1)+start <= isostartpos(ix)) && (cur_vertex(2)+start >= isostartpos(ix))
	to_draw = 1;
	%%%%%% DRAW EXONS
	if strand == '+' || ~param.flip
	  patch([isostartpos(ix) isostartpos(ix),...
		 cur_vertex(2)+start cur_vertex(2)+start],...
		[cnt-width, cnt+width, cnt+width, cnt-width],param.colors{1,5} );
	else
	  patch(-[isostartpos(ix) isostartpos(ix),...
		 cur_vertex(2)+start cur_vertex(2)+start],...
		[cnt-width, cnt+width, cnt+width, cnt-width],param.colors{1,5} );
	end
	%%%%%% ADD EXON LENGTH
	if strand == '+' || ~param.flip
	  h = text(mean([cur_vertex(2)+start isostartpos(ix)]), cnt-2*width, ...
		   sprintf('%i', cur_vertex(2)+start-isostartpos(ix)+1));
	else
	  h = text(-mean([cur_vertex(2)+start isostartpos(ix)]), cnt-2*width, ...
		   sprintf('%i', cur_vertex(2)+start-isostartpos(ix)+1));
	end
	set(h, 'HorizontalAlignment', 'center') ;
	set(h, 'FontSize',param.FontSize) ;
	set(h, 'FontWeight',param.FontWeight) ;
      end
    end
  end
end





  
function  cnt = viewgraph(vertices,edges,chr,strand,param,cnt)
%%%%%% ADD SUMMARY SPLICE GRAPH  
start = param.start;
stop = param.stop;
width = param.width;
width2 = param.width2;
step = param.step;

%%%%%% DRAW A LINE WITH DOTS WHERE N'S ARE.
cnt = cnt + 2*step;
hold on;
if strand == '+' || ~param.flip
  plot([param.start, param.stop],[cnt, cnt],'b');
else
  plot(-[param.start, param.stop],[cnt, cnt],'b');
end
str = load_mrna(chr,'+',param.start,param.stop);
nidx = find(lower(str)=='n');
if strand == '+' || ~param.flip
  plot(nidx+param.start,repmat(cnt,1,length(nidx)),'r.');
else
  plot(-(nidx+param.start),repmat(cnt,1,length(nidx)),'r.');
end
%if length(nidx) > 0, keyboard;end;

%%%%%% DRAW EXONS 
if ~isempty(vertices)
  cnt = cnt+1*step;
end
exon_num = zeros(1,stop-start);
exon_loc = zeros(1,stop-start);
exon_level = zeros(size(vertices,2));
for i=1:size(vertices,2)
  cur_vertex = vertices(:,i) - [start;start];
  exon_num(cur_vertex(1):cur_vertex(2)) = ...
      exon_num(cur_vertex(1):cur_vertex(2)) + 1;
  if all(exon_num<2)
    exon_loc = exon_num;
    level = 1;
  elseif max(exon_num)>size(exon_loc,1)
    exon_loc(end+1,:) = 0; 
    exon_loc(end,cur_vertex(1):cur_vertex(2)) = 1; 
    level = size(exon_loc,1);
  elseif max(exon_num)<=size(exon_loc,1)
    idx = min(find(all(exon_loc(:,cur_vertex(1):cur_vertex(2))==0,2))); 
    exon_loc(idx,cur_vertex(1):cur_vertex(2)) = 1;
    level = idx;
  end
 
  exon_level(i) = level;
  
  if strand == '+' || ~param.flip
    patch([cur_vertex(1)+start cur_vertex(1)+start,...
	   cur_vertex(2)+start cur_vertex(2)+start],...
	  [cnt+(level*step)-width, cnt+(level*step)+width,...
	   cnt+(level*step)+width, cnt+(level*step)-width],param.colors{1,3} );
  else
    patch(-[cur_vertex(1)+start cur_vertex(1)+start,...
	   cur_vertex(2)+start cur_vertex(2)+start],...
	  [cnt+(level*step)-width, cnt+(level*step)+width,...
	   cnt+(level*step)+width, cnt+(level*step)-width],param.colors{1,3} );
  end  
  if strand == '+' || ~param.flip
    h=text(mean([vertices(1,i),vertices(2,i)]), cnt+(level*step)-2*width, ...
	   sprintf('%i', abs(cur_vertex(2)-cur_vertex(1)+1))) ;
  else
    h=text(-mean([vertices(1,i),vertices(2,i)]), cnt+(level*step)-2*width, ...
	   sprintf('%i', abs(cur_vertex(2)-cur_vertex(1)+1))) ;
  end
  set(h, 'HorizontalAlignment', 'center') ;
  set(h, 'FontSize',param.FontSize) ;
  set(h, 'FontWeight',param.FontWeight) ;
end

%%%%%% DRAW INTRONS 
intron_loc = zeros(1,stop-start);
if (size(edges,1) > 1)

  for i=1:size(vertices,2)
    for j=i+1:size(vertices,2)
      if edges(i,j)
	if (vertices(1,i) < vertices(1,j))
	  istart = vertices(2,i);
	  istop =  vertices(1,j);
	  level1 = exon_level(i);
	  level2 = exon_level(j);
	else
	  istart = vertices(2,j);
	  istop =  vertices(1,i);
	  level1 = exon_level(j);
	  level2 = exon_level(i);
	end

	cur_intron = [istart-start,istop-start];
	intron_loc(cur_intron(1):cur_intron(2)) = ...
	    intron_loc(cur_intron(1):cur_intron(2)) + 1;
        leveli = (level1+level2).*.6;
	if strand == '+' || ~param.flip
	  h = line([istart, mean([istart,istop]),istop],...
		   [cnt+(level1*step)+width/2, cnt+((leveli+1)*step*0.75), cnt+(level2*step)+width/2]);
	else
	  h = line(-[istart, mean([istart,istop]),istop],...
		   [cnt+(level1*step)+width/2, cnt+((leveli+1)*step*0.75), cnt+(level2*step)+width/2]);
	end
        set(h, 'LineWidth',param.LineWidth) ;
        set(h, 'Color',param.colors{1,4}) ;

	if strand == '+' || ~param.flip
	  h = text(mean([istart istop]), cnt+((leveli+1)*step*0.75), sprintf('%i', abs(istop-istart)));
	else
	  h = text(-mean([istart istop]), cnt+((leveli+1)*step*0.75), sprintf('%i', abs(istop-istart)));
	end
	set(h, 'HorizontalAlignment', 'center') ;
	set(h, 'FontSize',param.FontSize) ;
	set(h, 'FontWeight',param.FontWeight) ;

      end
    end
  end
end
if ~isempty(exon_level)
  cnt = cnt+max(max(exon_level))*step;
end






function  cnt = viewname(name,chr,strand,param,cnt)
start = param.start;
stop = param.stop;
width = param.width;
width2 = param.width2;
step = param.step;


%%%%% DRAW GRID
if param.grid_on 
  gridsize = round((stop-start)/5);
  for i=ceil(start/gridsize)*gridsize:gridsize:floor(stop/gridsize)*gridsize,
    if strand == '+' || ~param.flip
      plot([i i], [step/4 cnt+1*step], 'k:') ;
      % plot([i i], [step/4 step*length(truepos)+step/2], 'k:') ;
    else
      plot(-[i i], [step/4 cnt+1*step], 'k:') ;
    end
    hold on
    if strand == '+' || ~param.flip
      h=text(i, step/8, sprintf('%i', i)) ;
    else
      h=text(-i, step/8, sprintf('%i', i)) ;
    end
    set(h, 'HorizontalAlignment', 'center') ;
    set(h, 'FontSize',param.FontSize) ;
    set(h, 'FontWeight',param.FontWeight) ;
  end ;
end
  
%%%%%% ADD TITLE 

if strand == '+' || ~param.flip
  %h = text(mean(start,stop),cnt+1.5*step, sprintf('%s\nChromosome %s %c',name,chr,strand)) ;
  h = text(mean(start,stop),cnt+1.5*step, sprintf('%s %c',strrep(chr,'_','\_'),strand)) ;
else
  %h = text(mean(start,stop),cnt+1.5*step, sprintf('%s\nChromosome %s %c',name,chr,strand)) ;
  h = text(-mean(start,stop),cnt+1.5*step, sprintf('%s %c',strrep(chr,'_','\_'),strand)) ;
end
set(h, 'HorizontalAlignment', 'left') ;
set(h, 'FontSize',param.FontSize ) ;
set(h, 'FontWeight',param.FontWeight) ;


if strand == '+' || ~param.flip
  axis([start-300 stop 0 cnt+2*step]) ;
else
  axis([-stop-300 -start 0 cnt+2*step]) ;
end
hold off
% axis([start stop 0 cnt+step*max((max([exon_loc,intron_loc])+0.5),3)]) ;
%keyboard

return

