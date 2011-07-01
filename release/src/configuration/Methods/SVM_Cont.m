function method = SVM_Cont(Content)

  
method = SVM(Content);
method.feature_type = 'char';
method.center = [];
method.kernels{1}.name = 'LINEAR';
method.genomewide_predict_fct = 'svm_content_predict';

for k=1:length(Content.wordlen)
  method.kernels{k}.wordlen = Content.wordlen{k};
  method.kernels{k}.stepping = Content.stepping ;
  method.kernels{k}.offset = Content.offset ;
end


