function include_preds(data)
% function include_preds(info,origpath,max_fd_rate)

info = data.ginfo;
origpath = data.origpath;
max_fd_rate = data.fd_rate;

%info = '~/altsplicedata/test/genome.config';origpath='/fml/ag-raetsch/home/ong/splicing';

addpath(sprintf('%s/utils', origpath)) ;
addpath ~/svn/tools/genomes
addpath ~/svn/tools/utils

% --- initialise genome
datadir = 'prediction';

genome_info = init_genome(info) ;
fprintf('Loading genome at: %s\n', genome_info.basedir) ;
origfilename=sprintf('%s/confirmed_sequences.mat', genome_info.basedir) ;
graphfilename=sprintf('%s/%s/cand_splicegraph_%02d.mat', genome_info.basedir, datadir, max_fd_rate) ;
if exist(graphfilename,'file')
  return
end

genes = load_cell(origfilename,'genes') ;
fprintf('loaded %s.\n', origfilename) ;


for ix = 1:length(genes)
  genes(ix).alt_exon = {};
  genes(ix).alt_intron = {};
  genes(ix).alt_5prime = {};
  genes(ix).alt_3prime = {};
end  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exon skip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf(1,'including predicted exon skips\n');
% exon skip, include new intron
eventfilename = sprintf('%s/%s/cand_alt_exon.mat', genome_info.basedir, datadir) ;
load(eventfilename,'cand_as_event');

for ix = 1:length(cand_as_event)
  if cand_as_event(ix).fd_rate > max_fd_rate/100
    continue;
  end
  gene = genes(cand_as_event(ix).genes_idx);
  if length(gene.splicegraph) == 2
    gene.splicegraph{3} = zeros(size(gene.splicegraph{1}));
    gene.splicegraph{4} = zeros(size(gene.splicegraph{2}));
  end
  exon1 = find(all(gene.splicegraph{1} == repmat(cand_as_event(ix).exon_5prime',1,size(gene.splicegraph{1},2))));
  exon2 = find(all(gene.splicegraph{1} == repmat(cand_as_event(ix).exon_3prime',1,size(gene.splicegraph{1},2))));
  assert(all((size(exon1)==[1,1]) & (size(exon2)==[1,1])))
  if gene.splicegraph{2}(exon1,exon2) == 1,
    fprintf('exon skip known\n');
    continue;
  end
  %exon1 = cand_as_event(ix).exon_idx(1);
  %exon2 = cand_as_event(ix).exon_idx(2);
  gene.splicegraph{2}(exon1,exon2) = 1;
  gene.splicegraph{2}(exon2,exon1) = 1;
  gene.splicegraph{4}(exon1,exon2) = cand_as_event(ix).fd_rate;
  gene.splicegraph{4}(exon2,exon1) = cand_as_event(ix).fd_rate;
  gene.alt_exon(end+1) = {[cand_as_event(ix).exon_idx]};
  
  genes(cand_as_event(ix).genes_idx) = gene;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intron retention
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'including predicted intron retention\n');
% intron retention, include new exon and join to old exons
eventfilename = sprintf('%s/%s/cand_alt_intron.mat', genome_info.basedir, datadir) ;
load(eventfilename,'cand_as_event');

for ix = 1:length(cand_as_event)
  if cand_as_event(ix).fd_rate > max_fd_rate/100
    continue;
  end
  gene = genes(cand_as_event(ix).genes_idx);
  if length(gene.splicegraph) == 2
    gene.splicegraph{3} = zeros(size(gene.splicegraph{1}));
    gene.splicegraph{4} = zeros(size(gene.splicegraph{2}));
  end
  if cand_as_event(ix).strand == '+'
    exon1 = find(all(gene.splicegraph{1} == repmat(cand_as_event(ix).exon_5prime',1,size(gene.splicegraph{1},2))));
    exon2 = find(all(gene.splicegraph{1} == repmat(cand_as_event(ix).exon_3prime',1,size(gene.splicegraph{1},2))));
    new_exon = [cand_as_event(ix).exon_5prime(1); cand_as_event(ix).exon_3prime(2)];
  else
    exon1 = find(all(gene.splicegraph{1} == repmat(cand_as_event(ix).exon_3prime',1,size(gene.splicegraph{1},2))));
    exon2 = find(all(gene.splicegraph{1} == repmat(cand_as_event(ix).exon_5prime',1,size(gene.splicegraph{1},2))));
    new_exon = [cand_as_event(ix).exon_3prime(1); cand_as_event(ix).exon_5prime(2)];
  end
  if ~(all((size(exon1)==[1,1]) & (size(exon2)==[1,1])))
    fprintf('multiple matches.\n');
    exon1 = exon1(1);
    exon2 = exon2(1);
  end
  if (cand_as_event(ix).exon_idx(3) ~= 0), 
    fprintf('known event.\n');
    assert(cand_as_event(ix).alt1_confirm_by < 0.2 * cand_as_event(ix).alt0_confirm_by);
    continue;
  end
  %assert(cand_as_event(ix).exon_idx(3) == 0);
  
  % grow the matrices by including the new exon after exon1
  gene.splicegraph{1} = [gene.splicegraph{1}(:,1:exon1),new_exon,gene.splicegraph{1}(:,exon1+1:end)];
  gene.splicegraph{3} = [gene.splicegraph{3}(:,1:exon1),repmat(cand_as_event(ix).fd_rate,2,1),gene.splicegraph{3}(:,exon1+1:end)];

  gene.splicegraph{2} = [gene.splicegraph{2}(1:exon1,1:exon1),gene.splicegraph{2}(1:exon1,exon1),...
		    gene.splicegraph{2}(1:exon1,exon1+1:end);...
		    gene.splicegraph{2}(exon1,1:exon1),zeros(1,exon2-exon1),...
		    gene.splicegraph{2}(exon2,exon2:end);...
		    gene.splicegraph{2}(exon1+1:end,1:exon1),...
		    [zeros(exon2-exon1-1,1);gene.splicegraph{2}(exon2:end,exon2)],...
		    gene.splicegraph{2}(exon1+1:end,exon1+1:end)];
  gene.splicegraph{4} = [gene.splicegraph{4}(1:exon1,1:exon1),gene.splicegraph{4}(1:exon1,exon1),...
		    gene.splicegraph{4}(1:exon1,exon1+1:end);...
		    gene.splicegraph{4}(exon1,1:exon1),zeros(1,exon2-exon1),...
		    gene.splicegraph{4}(exon2,exon2:end);...
		    gene.splicegraph{4}(exon1+1:end,1:exon1),...
		    [zeros(exon2-exon1-1,1);gene.splicegraph{4}(exon2:end,exon2)],...
		    gene.splicegraph{4}(exon1+1:end,exon1+1:end)];
  
  
  gene = shift_known_events(gene, exon1);
  gene.alt_intron(end+1) = {[exon1, exon2+1, exon1+1]};

  genes(cand_as_event(ix).genes_idx) = gene;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alternative 3 prime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'including predicted alternative 3prime\n');
% alt 3 prime, include new exon and join to old exon
eventfilename = sprintf('%s/%s/cand_alt_3prime.mat', genome_info.basedir, datadir) ;
load(eventfilename,'cand_as_event');


for ix = 1:length(cand_as_event)
  %if all(cand_as_event(ix).exon_3prime == cand_as_event(ix).exon_alternative)
  %  warning('alt_3prime: exon_3prime repeated.')
  %  continue;
  %end
  
  if cand_as_event(ix).fd_rate > max_fd_rate/100
    continue;
  end
  gene = genes(cand_as_event(ix).genes_idx);
  if length(gene.splicegraph) == 2
    gene.splicegraph{3} = zeros(size(gene.splicegraph{1}));
    gene.splicegraph{4} = zeros(size(gene.splicegraph{2}));
  end
  exon1 = find(all(gene.splicegraph{1} == repmat(cand_as_event(ix).exon_5prime',1,size(gene.splicegraph{1},2))));
  
  exon2a = find(all(gene.splicegraph{1} == repmat(cand_as_event(ix).exon_3prime',1,size(gene.splicegraph{1},2))));
  exon2b = find(all(gene.splicegraph{1} == repmat(cand_as_event(ix).exon_alternative',1,size(gene.splicegraph{1},2))));
  if ~(xor(all(size(exon2a)==[1,1]),all(size(exon2b)==[1,1]))) || ~(xor(isempty(exon2a),isempty(exon2b)))
    fprintf('alt_3prime: event %d in gene %d: ',ix,cand_as_event(ix).genes_idx);
    if ~isempty(exon2a) && ~isempty(exon2b)
      fprintf('known exons.\n');
      gene.splicegraph{2}(exon1,exon2b) = 1;
      gene.splicegraph{2}(exon1,exon2a) = 1;
      gene.splicegraph{4}(exon1,exon2b) = cand_as_event(ix).fd_rate;
      gene.splicegraph{4}(exon1,exon2a) = cand_as_event(ix).fd_rate;
      continue;
    else
      fprintf('multiple match.\n');
      if isempty(exon2a)
	exon2b = exon2b(1);
      end
      if isempty(exon2b)
	exon2a = exon2a(1);
      end
    end
  end
  
  assert(xor(all(size(exon2a)==[1,1]),all(size(exon2b)==[1,1])))
  assert(xor(isempty(exon2a),isempty(exon2b)))
  if isempty(exon2a)
    exon2 = exon2b;
    exon_new_pos = cand_as_event(ix).exon_3prime';
  else %isempty(exon2b)
    exon2 = exon2a;
    exon_new_pos = cand_as_event(ix).exon_alternative';
  end
  if ~(all((size(exon1)==[1,1]) & (size(exon2)==[1,1])))
    fprintf('multiple matches\n');
    exon1 = exon1(1);
    exon2 = exon2(1);
  end
  %exon1 = cand_as_event(ix).exon_idx(1);
  %exon2 = cand_as_event(ix).exon_idx(2);
  if exon2 < exon1
    temp = exon2;
    exon2 = exon1;
    exon1 = temp;
  end

  % put exon_alternative on the end and re-sort the exons.
  gene.splicegraph{1} = [gene.splicegraph{1},exon_new_pos];
  gene.splicegraph{3} = [gene.splicegraph{3},repmat(cand_as_event(ix).fd_rate,2,1)];
  

  if cand_as_event(ix).strand == '+'
    gene.splicegraph{2} = [gene.splicegraph{2},...
		    [zeros(exon1-1,1);1;zeros(exon2-exon1-1,1);gene.splicegraph{2}(exon2:end,exon2)];...
		    [zeros(1,exon1-1),1,zeros(1,exon2-exon1-1),gene.splicegraph{2}(exon2,exon2:end),0]];
    
    gene.splicegraph{4} = [gene.splicegraph{4},...
		    [zeros(exon1-1,1);cand_as_event(ix).fd_rate;zeros(exon2-exon1-1,1);gene.splicegraph{4}(exon2:end,exon2)];...
		    [zeros(1,exon1-1),cand_as_event(ix).fd_rate,zeros(1,exon2-exon1-1),gene.splicegraph{4}(exon2,exon2:end),0]];
  else
    m = size(gene.splicegraph{2},1);
    gene.splicegraph{2} = [gene.splicegraph{2},...
		    [gene.splicegraph{2}(1:exon1,exon1);zeros(exon2-exon1-1,1);1;zeros(m-exon2,1)];...
		    [gene.splicegraph{2}(exon1,1:exon1),zeros(1,exon2-exon1-1),1,zeros(1,m-exon2),0]];
    
    gene.splicegraph{4} = [gene.splicegraph{4},...
		    [gene.splicegraph{4}(1:exon1,exon1);zeros(exon2-exon1-1,1);cand_as_event(ix).fd_rate;zeros(m-exon2,1)];...
		    [gene.splicegraph{4}(exon1,1:exon1),zeros(1,exon2-exon1-1),cand_as_event(ix).fd_rate,zeros(1,m-exon2),0]];
  end
  
  [dummy,exon_order] = sort(gene.splicegraph{1}(1,:),2,'ascend');
  gene.splicegraph{1} = gene.splicegraph{1}(:,exon_order);
  gene.splicegraph{3} = gene.splicegraph{3}(:,exon_order);
  gene.splicegraph{2} = gene.splicegraph{2}(exon_order,exon_order);
  gene.splicegraph{4} = gene.splicegraph{4}(exon_order,exon_order);

  gene = shift_known_events(gene, find(exon_order==size(gene.splicegraph{1},2)));
  gene.alt_3prime(end+1) = {[find(exon_order==exon1), find(exon_order==exon2), ...
		    find(exon_order==size(gene.splicegraph{1},2))]};
  genes(cand_as_event(ix).genes_idx) = gene;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alternative 5 prime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'including predicted alternative 5prime\n');
% alt 5 prime, include new exon and join to old exon
eventfilename = sprintf('%s/%s/cand_alt_5prime.mat', genome_info.basedir, datadir) ;
load(eventfilename,'cand_as_event');


for ix = 1:length(cand_as_event)
  %if all(cand_as_event(ix).exon_5prime == cand_as_event(ix).exon_alternative)
  %  warning('alt_5prime: exon_5prime repeated.')
  %  continue;
  %end

  if cand_as_event(ix).fd_rate > max_fd_rate/100
    continue;
  end
  gene = genes(cand_as_event(ix).genes_idx);
  if length(gene.splicegraph) == 2
    gene.splicegraph{3} = zeros(size(gene.splicegraph{1}));
    gene.splicegraph{4} = zeros(size(gene.splicegraph{2}));
  end
  
  exon2 = find(all(gene.splicegraph{1} == repmat(cand_as_event(ix).exon_3prime',1,size(gene.splicegraph{1},2))));
  exon1a = find(all(gene.splicegraph{1} == repmat(cand_as_event(ix).exon_5prime',1,size(gene.splicegraph{1},2))));
  exon1b = find(all(gene.splicegraph{1} == repmat(cand_as_event(ix).exon_alternative',1,size(gene.splicegraph{1},2))));
  if ~(xor(all(size(exon1a)==[1,1]),all(size(exon1b)==[1,1]))) || ~(xor(isempty(exon1a),isempty(exon1b)))
    fprintf('alt_5prime: event %d in gene %d: ',ix,cand_as_event(ix).genes_idx);
    if ~isempty(exon1a) && ~isempty(exon1b)
      fprintf('known exons\n');
      gene.splicegraph{2}(exon2,exon1b) = 1;
      gene.splicegraph{2}(exon2,exon1a) = 1;
      gene.splicegraph{4}(exon2,exon1b) = cand_as_event(ix).fd_rate;
      gene.splicegraph{4}(exon2,exon1a) = cand_as_event(ix).fd_rate;
      continue;
    else
      fprintf('multiple match.\n');
      if isempty(exon1a)
	exon1b = exon1b(1);
      end
      if isempty(exon1b)
	exon1a = exon1a(1);
      end
    end
  end
  
  assert(xor(all(size(exon1a)==[1,1]),all(size(exon1b)==[1,1])))
  assert(xor(isempty(exon1a),isempty(exon1b)))
  if isempty(exon1a)
    exon1 = exon1b;
    exon_new_pos = cand_as_event(ix).exon_5prime';
  else %isempty(exon1b)
    exon1 = exon1a;
    exon_new_pos = cand_as_event(ix).exon_alternative';
  end    
  exon2 = find(all(gene.splicegraph{1} == repmat(cand_as_event(ix).exon_3prime',1,size(gene.splicegraph{1},2))));
  if ~(all((size(exon1)==[1,1]) & (size(exon2)==[1,1])))
    fprintf('multiple matches\n');
    exon1 = exon1(1);
    exon2 = exon2(1);
  end
  %exon1 = cand_as_event(ix).exon_idx(1);
  %exon2 = cand_as_event(ix).exon_idx(2);
  if exon2 < exon1
    temp = exon2;
    exon2 = exon1;
    exon1 = temp;
  end

  % put exon_alternative on the end and re-sort the exons.
  gene.splicegraph{1} = [gene.splicegraph{1},exon_new_pos];
  gene.splicegraph{3} = [gene.splicegraph{3},repmat(cand_as_event(ix).fd_rate,2,1)];
  

  if cand_as_event(ix).strand == '+'
    m = size(gene.splicegraph{2},1);
    gene.splicegraph{2} = [gene.splicegraph{2},...
		    [gene.splicegraph{2}(1:exon1,exon1);zeros(exon2-exon1-1,1);1;zeros(m-exon2,1)];...
		    [gene.splicegraph{2}(exon1,1:exon1),zeros(1,exon2-exon1-1),1,zeros(1,m-exon2),0]];
    
    gene.splicegraph{4} = [gene.splicegraph{4},...
		    [gene.splicegraph{4}(1:exon1,exon1);zeros(exon2-exon1-1,1);cand_as_event(ix).fd_rate;zeros(m-exon2,1)];...
		    [gene.splicegraph{4}(exon1,1:exon1),zeros(1,exon2-exon1-1),cand_as_event(ix).fd_rate,zeros(1,m-exon2),0]];
  else
    gene.splicegraph{2} = [gene.splicegraph{2},...
		    [zeros(exon1-1,1);1;zeros(exon2-exon1-1,1);gene.splicegraph{2}(exon2:end,exon2)];...
		    [zeros(1,exon1-1),1,zeros(1,exon2-exon1-1),gene.splicegraph{2}(exon2,exon2:end),0]];
    
    gene.splicegraph{4} = [gene.splicegraph{4},...
		    [zeros(exon1-1,1);cand_as_event(ix).fd_rate;zeros(exon2-exon1-1,1);gene.splicegraph{4}(exon2:end,exon2)];...
		    [zeros(1,exon1-1),cand_as_event(ix).fd_rate,zeros(1,exon2-exon1-1),gene.splicegraph{4}(exon2,exon2:end),0]];
  end
  
  [dummy,exon_order] = sort(gene.splicegraph{1}(1,:),2,'ascend');
  gene.splicegraph{1} = gene.splicegraph{1}(:,exon_order);
  gene.splicegraph{3} = gene.splicegraph{3}(:,exon_order);
  gene.splicegraph{2} = gene.splicegraph{2}(exon_order,exon_order);
  gene.splicegraph{4} = gene.splicegraph{4}(exon_order,exon_order);

  gene = shift_known_events(gene, find(exon_order==size(gene.splicegraph{1},2)));
  gene.alt_5prime(end+1) = {[find(exon_order==exon1), find(exon_order==exon2), ...
		    find(exon_order==size(gene.splicegraph{1},2))]};
  genes(cand_as_event(ix).genes_idx) = gene;  
end


fprintf(1,'Saving %s...',graphfilename);
save_cell(graphfilename, genes, 'genes');
fprintf(1, 'done\n');












return





function gene = shift_known_events(gene, exonpos)

for ixa = 1:length(gene.alt_exon)
  incridx = find(gene.alt_exon{ixa}>exonpos);
  newevent = gene.alt_exon{ixa};
  newevent(incridx) = newevent(incridx) + 1;
  gene.alt_exon(ixa) = {newevent};
end
for ixa = 1:length(gene.alt_intron)
  incridx = find(gene.alt_intron{ixa}>exonpos);
  newevent = gene.alt_intron{ixa};
  newevent(incridx) = newevent(incridx) + 1;
  gene.alt_intron(ixa) = {newevent};
end
for ixa = 1:length(gene.alt_5prime)
  incridx = find(gene.alt_5prime{ixa}>exonpos);
  newevent = gene.alt_5prime{ixa};
  newevent(incridx) = newevent(incridx) + 1;
  gene.alt_5prime(ixa) = {newevent};
end
for ixa = 1:length(gene.alt_3prime)
  incridx = find(gene.alt_3prime{ixa}>exonpos);
  newevent = gene.alt_3prime{ixa};
  newevent(incridx) = newevent(incridx) + 1;
  gene.alt_3prime(ixa) = {newevent};
end

return







% some code to view the splice graphs with predicted events
if false
  load ~/altsplicedata/test/prediction/cand_splicegraph.mat
  for ix = 1:length(genes)
    if length(genes(ix).splicegraph) == 4
      fprintf('gene id:%d; strand=%c\n',ix,genes(ix).strands(1));
      fprintf('alt_exon\n');
      for ixa = 1:length(genes(ix).alt_exon)
	genes(ix).alt_exon{ixa}
      end
      fprintf('alt_intron\n');
      for ixa = 1:length(genes(ix).alt_intron)
	genes(ix).alt_intron{ixa}
      end
      fprintf('alt_5prime\n');
      for ixa = 1:length(genes(ix).alt_5prime)
	genes(ix).alt_5prime{ixa}
      end
      fprintf('alt_3prime\n');
      for ixa = 1:length(genes(ix).alt_3prime)
	genes(ix).alt_3prime{ixa}
      end

      viewsplicegraph_preds(genes(ix),1);
      pause;
    end
  end
end
