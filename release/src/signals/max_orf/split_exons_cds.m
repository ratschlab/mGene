function [cds_exons, utr5_exons, utr3_exons, ok] = split_exons_cds(exons, tis_pos, cdsStop_pos, strand) ;
% [cds_exons, utr5_exons, utr3_exons, ok] = split_exons_cds(exons, tis_pos, cdsStop_pos, strand) ;

cds_exons = [] ;
utr5_exons = [] ;
utr3_exons = [] ;
ok=1 ;

cds_exons = exons ;
tis_exon = find(cds_exons(:,1)<=tis_pos & cds_exons(:,2)>=tis_pos) ;
if length(tis_exon)~=1,
  ok=0 ;
  return ;
end ;

if strand == '+',
  utr5_exons = cds_exons(1:tis_exon,:) ;
  utr5_exons(end,2) = tis_pos ;
  if utr5_exons(end,2)-utr5_exons(end,1)==0,
    utr5_exons(end,:) = [] ;
  elseif utr5_exons(end,2)-utr5_exons(end,1)<0,
    assert(0) % this should not happen
    keyboard
    utr5_exons(end,:) = [] ;
    ok = 0 ;
    return ;
  end ;
  
  cds_exons = cds_exons(tis_exon:end, :) ;
  cds_exons(1,1) = tis_pos ;
else
  utr5_exons = cds_exons(tis_exon:end, :) ;
  utr5_exons(1,1) = tis_pos ;
  if utr5_exons(1,2)-utr5_exons(1,1)==0,
    utr5_exons(1,:) = [] ;
  elseif utr5_exons(1,2)-utr5_exons(1,1)<0,
    assert(0) % this should not happen
    keyboard
    utr5_exons(1,:) = [] ;
    ok = 0 ;
    return ;
  end ;
  
  cds_exons = cds_exons(1:tis_exon, :) ;
  cds_exons(end,2) = tis_pos ;
end ;

stop_exon = find(cds_exons(:,1)<=cdsStop_pos & cds_exons(:,2)>=cdsStop_pos) ;
if length(stop_exon)~=1,
  %assert(0) % this should not happen
  %keyboard
  ok=0 ;
  return ;
end ;

if strand == '+',
  utr3_exons = cds_exons(stop_exon:end, :) ;
  utr3_exons(1,1) = cdsStop_pos ;
  if utr3_exons(1,2)-utr3_exons(1,1)==0,
    utr3_exons(1,:) = [] ;
  elseif utr3_exons(1,2)-utr3_exons(1,1)<0,
    assert(0) % this should not happen
    utr3_exons(1,:) = [] ;
    keyboard
  end ;
  
  cds_exons = cds_exons(1:stop_exon, :) ;
  cds_exons(end,2) = cdsStop_pos ;
else
  utr3_exons = cds_exons(1:stop_exon,:) ;
  utr3_exons(end,2) = cdsStop_pos ;
  if utr3_exons(end,2)-utr3_exons(end,1)==0,
    utr3_exons(end,:) = [] ;
  elseif utr3_exons(end,2)-utr3_exons(end,1)<0,
    assert(0) % this should not happen
    keyboard
    utr3_exons(end,:) = [] ;
  end ;
  
  cds_exons = cds_exons(stop_exon:end, :) ;
  cds_exons(1,1) = cdsStop_pos ;
end ;

ok=1 ;
