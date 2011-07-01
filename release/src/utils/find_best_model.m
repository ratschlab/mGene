function fn_train_best = find_best_model(output_dir)

max_res = 0;
max_iter = 0;
train_name = sprintf('%s/lsl/train/train_iteration%i.mat', output_dir, 1)
C = load_C(train_name);

for j = 1:300
  %fn_pred = sprintf('%s/lsl/predictions/iteration_%i/results.mat', output_dir, j);
  fn_train = sprintf('%s/lsl/train/train_iteration%i.mat', output_dir, j);
  %if fexist(fn_pred)
  if fexist(fn_train)
    train_name = sprintf('%s/lsl/train/train_iteration%i.mat', output_dir, j);
	max_iter = j;
    %l = load(fn_pred);
    %current_C = load_C(train_name);
    % if this is just a nother iteration with the same model 
    % then use the converged predictor even if the result on 
    % the test set is worse 
    if 0%l.res.cds_transcripts.F>max_res || current_C==C
      max_res = l.res.cds_transcripts.F;
      fn_pred_best = fn_pred;
      max_iter = j;
      C = current_C;
    end
  end 
end
fn_train_best = sprintf('%s/lsl/train/train_iteration%i.mat', output_dir, max_iter)
if ~(fexist(fn_train_best))
	error('file not found: %s', fn_train_best);
end

return

function C= load_C(train_file)

  load(train_file, 'QP')
  C = QP.Q(1,1);
return
