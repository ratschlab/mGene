function [blocks, QP, viterbi_nbest,num_new_constraints, time] = add_constraints(blocks, block_idx, num_examples, PAR, QP, viterbi_nbest);
%   [blocks, QP, viterbi_nbest] = add_constraints(blocks,block_idx,num_examples,PAR,QP,viterbi_nbest);
% 
%
if ~isfield(PAR, 'fn_log')
  %logfile = '~/tmp/logfile'
  fid = 1;%fopen(logfile, 'a+');
else
  fid = fopen(PAR.fn_log,'a+') ;
end

verb_level=1 ;

num_new_constraints = 0;
max_neg  = PAR.LSL.method.max_neg ;
opt_step = PAR.LSL.method.opt_step;
model = PAR.model ;
time = [];

clear A1 A1_sparse
A1 = zeros(opt_step,model.cnt_parameters) ;
A2 = sparse([], [], [], length(blocks)*max_neg, size(QP.A,2)-model.cnt_parameters, length(blocks)*max_neg) ; 
b = zeros(opt_step, 1) ;
q = 0 ; r = 0 ; q2 = 0 ;


blocks = fix_blocks(blocks, model) ;
model = fix_model(model) ;

[xis, parameters] = res2param(QP.res, num_examples, model, QP, 0, PAR.weight_names);
param_vec = weights2vector(parameters, PAR.weight_names) ;
num_correct = 0 ;
fprintf('\n') ;
losses = [] ;

w_fields = PAR.weight_names;

% check the weights from time to time
check_compute_weights=(rand(1)>0.95) ;
for id = 1:length(blocks),
  %blocks(id).pred_use_approximation = 0 ;
  block = blocks(id) ;

  cnt = 1; 
  prediction_fn = sprintf('%spred_block%i.mat', PAR.FN.output_lsl.fn_train_block_preds, block.id);
  while ~fexist(prediction_fn)
    fprintf('.') ;
    pause(10)
    if mod(cnt, 10)==0
      fprintf(fid,'\nwaiting for prediction: block: %i, example %i\n ',block.id,id);
    end
	if mod(cnt, 20)==0
      fprintf(fid,'\nignoring prediction: block: %i, example %i\n ',block.id,id);
	  break;
	end
	cnt = cnt + 1;
  end
  if (cnt>=10)
    continue;
  end


  load(prediction_fn, 'pred', 'decode_time', 'true_path_weights'); 
  assert(isfield(pred, 'loss')) ;

  if isempty(pred), warning('empty'); continue; end ;
  if exist('decode_time'), 
    time(id) = decode_time;
    clear decode_time;
  end
  all_identical = 1 ;
  for pred_num = 1:length(pred)
    pred(pred_num).weights = reorder_fields(pred(pred_num).weights, w_fields);
    
    if 0%check_compute_weights,%this check cant work if the approximation with reduced candidate list is used combined with the long_transitions
      [pred2(pred_num).weights, loss2(pred_num)] = compute_weights(pred(pred_num), block, model, PAR.FN, parameters) ;    
      fnw = fieldnames(pred2(pred_num).weights) ;
      for kk=1:length(fnw),
        diff = norm(pred2(pred_num).weights.(fnw{kk})-pred(pred_num).weights.(fnw{kk})) ;
        if diff>1e-4,
          fprintf('compute weights diff: %s => %1.2f\n', fnw{kk}, diff) ;
          if keyboard_allowed(),
            keyboard ;
          end ;
        end ;
      end ;
    end ;

    loss(pred_num) = pred(pred_num).loss ;

    identical = 0 ;
    for t=1:length(block.truth)
      identical = equal(pred(pred_num).path, block.truth(t).path) & equal(pred(pred_num).pos, block.truth(t).pos) ;
      identical = identical | (loss(pred_num)==0) ;
      if identical, break ; end ;
    end ;
    if ~identical, all_identical=0; end ;

    %vec = weights2vector(block.truth(1).weights, PAR.weight_names) - weights2vector(pred(pred_num).weights, PAR.weight_names);
    % if the approximation (with less candidate positions) is used in combination with long_transitions then the weights have to 
    % be recomputed by gen_paths. the reason is that the predictions with long_transitions differ a little if the 
    % candidate lists are not the same 
    vec = weights2vector(true_path_weights, PAR.weight_names) - weights2vector(pred(pred_num).weights, PAR.weight_names);
    diff_before = -sum(param_vec.*vec);
    
    if identical,
      if loss(pred_num)>0,
        warning('unexpected: loss>0') ;
        if keyboard_allowed(), keyboard ; end ;
        %assert(~loss(pred_num)) ;
      end ;
      if pred_num==1, 
        num_correct = num_correct+1;
      end
      if verb_level>2,
        fprintf(fid, 'disabling gen_path approximation\n') ;
      end ;
      block.pred_use_approximation = 0 ;
    else
      all_identical=0 ;
      assert(norm(vec)>1e-8) ;
      %if loss(pred_num==0) && keyboard_allowed, keyboard; end 
      assert(loss(pred_num)>0) ;
      
      if diff_before+loss(pred_num)-xis(block_idx(id)) < 1e-7
        %fprintf(fid,'no margin violation') ;
        %fprintf('\nno margin violation\n') ;
        if verb_level>1,
          fprintf('No margin violation: %f > %f = %f - %f\n', -diff_before, loss(pred_num) - xis(block_idx(id)), loss(pred_num), xis(block_idx(id))) ;
        end ;
        if -diff_before - loss(pred_num)>1e-3 && block.pred_use_approximation==0,
  	  warning('true path scores higher than max_path (%1.4f)', -diff_before - loss(pred_num))
   	  if keyboard_allowed,
            keyboard
          end 
        end
        if pred_num==length(pred)
          if verb_level>2,
            fprintf(fid,'no constraint added\n');
          end ;
        end
        if verb_level>2,
          fprintf(fid, 'disabling gen_path approximation\n') ;
        end ;
	block.pred_use_approximation = 0 ;
        continue
      end
      if verb_level>1,
        fprintf('Adding violated constraint %f > %f = %f - %f \n', -diff_before, loss(pred_num) - xis(block_idx(id)), loss(pred_num), xis(block_idx(id))) ;
      end ;
      num_new_constraints = num_new_constraints+1;
      q  = q+1 ; 
      q2 = q2+1 ;
      A1(q2,1:model.cnt_parameters) = -vec;
      b(q2)= -loss(pred_num) ;
      A2(q, block_idx(id)) = -1 ; % slack
      
      if q2==size(A1,1)
        fprintf('*\n') ;  
        A1 = sparse(A1(1:q2,:)) ;
        assert(q==size(A1,1));
        A2 = A2(1:q,:) ;
        A = [A1 A2] ;
        QP.A = [QP.A; A] ;
        QP.b = [QP.b; b] ;
        
        %tic;[res,lambda,how] = qp_solve(QP.lpenv, QP.Q, QP.f, QP.A, QP.b, QP.lb, QP.ub, 0, 1, PAR.solver);toc;
        %tic;[res,lambda,how] = qp_solve(QP.lpenv, QP.Q, QP.f, QP.A, QP.b, QP.lb, QP.ub, 0, 1, 'sift');toc;
        
        %OBJ = sum(res.*f) + 0.5*res'*QP.Q*res ;
        %[xis,parameters,hvar,H] = res2param(res,num_examples,model);
        %sum_xis = sum(xis) ;

        q2 = 0; q = 0 ;
        A1 = zeros(opt_step,model.cnt_parameters) ;
        A2 = sparse([],[],[],length(blocks)*max_neg, size(QP.A, 2)-model.cnt_parameters, length(blocks)*max_neg) ; 
        b = zeros(opt_step, 1) ;
      end ;
      if isempty(block.pred)
        block.pred = pred(pred_num) ; 
      else
        if ~isfield(block.pred, 'loss')
          block.pred(end).loss=nan ;
          block.pred=orderfields(block.pred) ;
        end ;
        if length(block.pred)==1,
          block.pred = orderfields(block.pred) ;
        end ;
        block.pred(end+1) = pred(pred_num) ; 
      end ;
      block.pred(end).loss = loss(pred_num) ;

      %fprintf(fid, 'enabling gen_path approximation\n') ;
      %block.pred_use_approximation = 1 ;
    end ;
  end ; %%% loop over predictions 
  %if all_identical,
    %viterbi_nbest(id)=viterbi_nbest(id)+1 ;
    %fprintf('increasing number of nbest path to %i\n', viterbi_nbest(id)) ;
  %end ;
  blocks(id) = block ;
  losses(id) = mean(loss) ;

end ;%%% loop over blocks

if verb_level>0,
  if fid~=1,
    fprintf(fid,'id=%i    correct=%i(%2.2f)  average_loss=%2.2f  num_new_constraints=%i\r', id,num_correct, num_correct/id, mean(losses), num_new_constraints);
  end ;
  fprintf('\nid=%i    correct=%i(%2.2f)  average_loss=%2.2f  num_new_constraints=%i\n', id, num_correct,num_correct/id, mean(losses), num_new_constraints);
end ;

A1     = sparse(A1(1:q2,:)) ;
b      = b(1:q2) ;
A2     = A2(1:q,:) ;
A      = [A1 A2] ;
QP.A   = [QP.A;A] ;
QP.b   = [QP.b;b] ;

if fid>1
  fclose(fid) ;
end

% it does not help anyways (GR)
viterbi_nbest(:)=1 ;

return

%%% debug stuff

%%% check if weights of true path are correct
bb = load_struct('/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/p_pacificus/lsl/trsk_nc/output/lsl/data/training_blocks.mat', 'blocks')
bidx = find([bb.id]==block.id);
block2 = bb(bidx);
[weights2,loss2,score2,losses2,scores2] = compute_weights(block2.truth(1), block2, model, PAR.FN,[],1);
isequal(block.truth.weights, weights2)

block.truth.segments(:, 2)-block.truth.segments(:, 1);


