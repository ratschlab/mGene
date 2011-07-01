function res = compute_SN_SP(res,signals,end_tol);
  
names = {'nucleotides','cds_nucleotides','exons','cds_exons','transcripts','cds_transcripts','genes'};  
for i=1:length(names)
  level = res.(names{i});
  level.SN = level.num_corr/level.num_obsv;
  level.SP = level.num_corr/level.num_pred;
  level.mean_SNSP = [mean([level.SN,level.SP])] ;
  % harmonic mean: F-measure
  level.F = 2*level.SN*level.SP/(level.SN+level.SP);
  if isfield(level,'ME')
    level.ME_gene_id =unique(level.ME_gene_id);
    level.WE_gene_id =unique(level.WE_gene_id);
  end
  res.(names{i}) = level;
end

names = {'tis','cdsStop','tss','cleave','acc','don'};

for i=1:length(names)
  level = signals.(names{i});
  level.SN = level.num_corr/level.num_obsv;
  level.SP = level.num_corr/level.num_pred;
  level.mean_SNSP = [mean([level.SN,level.SP])] ;
  % harmonic mean: F-measure
  level.F = 2*level.SN*level.SP/(level.SN+level.SP);
  res.(names{i}) = level;
end


tss = res.tss;
cleave = res.cleave;

for i=end_tol
  tss.(['SN' num2str(i)]) = tss.(['num_corr' num2str(i)])/tss.num_obsv;
  tss.(['SP' num2str(i)]) = tss.(['num_corr' num2str(i)])/tss.num_pred;
  tss.(['mean_SNSP' num2str(i)]) = mean([tss.(['SN' num2str(i)]),tss.(['SP' num2str(i)])]); 
  
  cleave.(['SN' num2str(i)]) = cleave.(['num_corr' num2str(i)])/cleave.num_obsv;
  cleave.(['SP' num2str(i)]) = cleave.(['num_corr' num2str(i)])/cleave.num_pred;
  cleave.(['mean_SNSP' num2str(i)]) = mean([cleave.(['SN' num2str(i)]),cleave.(['SP' num2str(i)])]); 
end
                      
res.tss = tss;
res.cleave = cleave;

