function [neg pos] = fit_length_hist(pos, neg, steps, reduce_pos, make_plot)

% pos and neg need to be cell arrays of strings

if nargin<5
  make_plots=0;
end

if nargin<4
  reduce_pos=0;
end

% count lengths 
neg_lengths = zeros(1,length(neg));
for j=1:length(neg)
  neg_lengths(j)=length(neg{j});
end
pos_lengths = zeros(1,length(pos));
for j=1:length(pos)
  pos_lengths(j)=length(pos{j});
end

max_len = max([pos_lengths neg_lengths]);
step=ceil(max_len/steps);
bins = logspace(0,log10(max_len),steps);
bins(1)=0;
if max(bins)<max_len,
  bins = [bins max_len];
end

pos_hist = myhistc(pos_lengths,bins);
neg_hist = myhistc(neg_lengths,bins);
if make_plots
  figure;
  plot(1:length(bins),neg_hist,1:length(bins),pos_hist)
end

% scale histograms such that they are comparable
neg_hist_norm = neg_hist/length(neg);
pos_hist_norm = pos_hist/length(pos);

if make_plots
  figure;
  plot(1:length(bins),neg_hist_norm,1:length(bins),pos_hist_norm)
end

% calculate the number of sequences that has to be deleted for each length bin
overh = ceil(length(neg)*(neg_hist_norm-pos_hist_norm));

rm_idx=[];
for j = 1:length(neg)  
  len = length(neg{j});
  bin = sum(bins<=len);
  if overh(bin)>0
    overh(bin) = overh(bin)-1;
    rm_idx = [rm_idx j];
  end
end
waste = {neg{rm_idx}};
w_lens = neg_lengths(rm_idx);
neg(rm_idx)=[];
neg_lengths_2 = neg_lengths;
neg_lengths_2(rm_idx)=[];

% sort waste
[w_lens idx] = sort(w_lens,'descend');
waste = {waste{idx}};
clear idx

% recycle
new_neg = {}; 
j = 0;
while j<length(waste)
  j = j+1;
  %idcs = find(overh<0); 
  %idx = idxs(ceil(rand*length(idcs)));
  [tmp idx] = min(overh);
  if tmp>=0
    break;
  elseif w_lens(j)<bins(idx)
    continue;
  elseif idx == length(bins)
    continue;
  else
    len = floor(bins(idx)+rand*(bins(idx+1)-bins(idx)));
    p = min(length(waste{j}),len);
    new_neg = {new_neg{:} waste{j}(1:p)};
    if p<w_lens(j)-30 
      waste = {waste{:} waste{j}(p+1:end)};
      w_lens(end+1) = length(waste{end});
    end
    overh(idx) = overh(idx)+1;
  end
end

neg_hist_2 = myhistc(neg_lengths_2,bins);

if make_plots
  figure;
  plot(1:length(bins),neg_hist_2,1:length(bins),pos_hist)
end

neg = {neg{:} new_neg{:}};

if make_plots
  for j=1:length(neg)
    neg_lengths_3(j)=length(neg{j});
  end
  neg_hist_3 = myhistc(neg_lengths_3,bins);

  hist_opt = neg_hist_3-overh;
  figure;
  plot(1:length(bins),neg_hist_3/sum(neg_hist_3),1:length(bins),pos_hist/sum(pos_hist))
end

if reduce_pos
  rm_idx=[];
  for j = 1:length(pos)  
    len = length(pos{j});
    bin = sum(bins<=len);
    if overh(bin)<0
      overh(bin) = overh(bin)+1;
      rm_idx = [rm_idx j];
    end
  end
  pos(rm_idx) = [];
end
%handle = gca;
%prop = get(handle);
%prop.XTick = bins;
%set(handle, 'XTick',bins);

%keyboard 



