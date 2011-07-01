function calc_plifs(PAR_1)
fid=1;
q=0;

fn_SVMs = PAR_1.FN.output_sig.(PAR_1.Signal_name).fn_SVMs;
method = PAR_1.Signals.(PAR_1.Signal_name).method;
num_partitions = PAR_1.SETs.num_partitions;

for p = 1:num_partitions
  fprintf(fid,'\n\n-----------\n');
  fprintf(fid,'PARTITION %i\n',p);
  fprintf(fid,'-----------\n');
  fn_SVM = sprintf('%sbest_partition=%i.mat', fn_SVMs,p) ;

  %--------------------------------------------------------------------------------%
  % load parameters from SVM file
  %--------------------------------------------------------------------------------%
  load(fn_SVM,'Train');
  PAR = Train.PAR;
  PAR.fn_SVMs=fn_SVM;
  clear Train
  
  %--------------------------------------------------------------------------------%
  % learn PLIFS on evaluation set and convert output on test set
  %--------------------------------------------------------------------------------%
  PAR.fn_pred = sprintf('%spartition=%i.mat',fn_SVMs,p);
  PAR.fn_plifs = sprintf('%spartition=%i_plifs.mat',fn_SVMs,p);
  assert(fexist(PAR.fn_pred)==1)
  assert(fexist(PAR.fn_plifs)==1)
  fprintf(fid,'\n already predicted on set ');
  fprintf(fid,'%i',PAR.SETs.which_pred);
        
  %--------------------------------------------------------------------------------%
  %%% learn PLIFS on evaluation set
  %--------------------------------------------------------------------------------%
  clear  prob prob_cum limits Conf Conf_cum
  fprintf('\nload fn_plif\n')
  load(PAR.fn_plifs,'LT', 'output');
  if ~all(ismember({'Conf','Conf_cum'},who_file(PAR.fn_plifs)))
    if ~all(ismember({'prob','prob_cum', 'limits'},who_file(PAR.fn_SVMs)))
      disp('learn PLIFs')
      assert(length(output)==length(LT))
      [prob, prob_cum,limits] = learn_plifs(output,LT,method.plif.bins,1,method.plif.min_incr );
      %save(PAR.fn_SVMs,'-append','prob', 'prob_cum','limits');
      save_append(PAR.fn_SVMs, 1, 'prob', prob, 'prob_cum', prob_cum, ...
                  'limits', limits);
    else
      load(PAR.fn_SVMs,'prob', 'prob_cum','limits');
    end
    disp('Convert SVM Output to Confidence on fn_plif')
    [Conf, Conf_cum] = Out2Conf(output,prob, prob_cum,limits);
    %save(PAR.fn_plifs,'-append','Conf','Conf_cum');
    save_append(PAR.fn_plifs,1, 'Conf', Conf, 'Conf_cum', Conf_cum);
    clear output Conf Conf_cum LT
  end
  %--------------------------------------------------------------------------------%
  %%% Convert SVM Outputs on test set 
  %--------------------------------------------------------------------------------%
  load(PAR.fn_pred,'output');
  if ~all(ismember({'Conf','Conf_cum'},who_file(PAR.fn_pred)))
	load(PAR.fn_SVMs,'prob', 'prob_cum','limits');
    disp('Convert SVM Output to Confidence on fn_pred')
    [Conf, Conf_cum] = Out2Conf(output,prob, prob_cum,limits);
    %save(PAR.fn_pred,'-append','Conf','Conf_cum');
    save_append(PAR.fn_pred, 1,'Conf', Conf, 'Conf_cum', Conf_cum);
  end
end%loop over partitions
% [jobinfo,num_crashed] = rproc_wait(jobinfo,20,1,0);
