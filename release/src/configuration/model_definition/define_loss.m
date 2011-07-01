function model = define_loss(model)
  

% create loss matrices

equiv = [];
  
different = [model.segments.intergenic, model.segments.cds_exon, 3 ;
             model.segments.intergenic, model.segments.intron, 3;
			 model.segments.intergenic, model.segments.utr5exon, 2;
			 model.segments.intergenic, model.segments.utr3exon, 2;
             model.segments.intron, model.segments.cds_exon, 2];
if isfield(model.use.segments, 'rna_seq_polya')&&model.use.segments.rna_seq_polya
	different = [different;
				model.segments.rna_seq_polya, model.segments.cds_exon, 3;
				model.segments.rna_seq_polya, model.segments.intergenic, 3;
				model.segments.rna_seq_polya, model.segments.utr3exon, 0];
end
%--------------------------------
%if ~model.use.loss.orf_only
%  different = [different;
%               model.segments.intergenic, model.segments.utr5exon, 2;
%               model.segments.intergenic, model.segments.utr3exon, 2];
%  if model.use.segments.intercistronic
%    different = [different;
%                 model.segments.intercistronic, model.segments.utr5exon, 2;
%                 model.segments.intercistronic, model.segments.utr3exon, 2;
%                 model.segments.intergenictrans, model.segments.utr5exon, 2;
%                 model.segments.intergenictrans, model.segments.utr3exon, 2];
%    
%  end
%else
%  equiv = [equiv;
%           model.segments.intergenic model.segments.utr3exon ;
%           model.segments.intergenic model.segments.utr5exon] ;
%  if model.use.segments.intercistronic
%    equiv = [equiv;
%             model.segments.intercistronic model.segments.utr3exon ; 
%             model.segments.intercistronic model.segments.utr5exon ;
%             model.segments.intergenictrans model.segments.utr3exon ;
%             model.segments.intergenictrans model.segments.utr5exon] ;
%  end
%end
%%--------------------------------
%if model.use.segments.polya_tail
%  equiv = [equiv;
%           model.segments.utr3exon, model.segments.polya_tail];
%  different = [different;
%               model.segments.intergenic, model.segments.polya_tail, 1];    
%end
% %-------------------------------- 
%if model.use.segments.intercistronic
%  
%  different = [different;
%               model.segments.intercistronic, model.segments.cds_exon, 3 ;
%               model.segments.intercistronic, model.segments.intron, 3;
%               model.segments.intercistronic, model.segments.transexon, 2;
%               model.segments.intergenic, model.segments.intercistronic, 0.1;
%               model.segments.intergenictrans, model.segments.intercistronic, 0.1;
%               model.segments.intergenictrans, model.segments.cds_exon, 3 ;
%               model.segments.intergenictrans, model.segments.intron, 3;
%               model.segments.intergenictrans, model.segments.transexon, 2;
%               model.segments.intergenic, model.segments.intergenictrans, 0.1];
%
%  if model.use.segments.polya_tail
%    different = [different;
%                 model.segments.intercistronic, model.segments.polya_tail,1;
%                 model.segments.intergenictrans, model.segments.polya_tail, 1];
%  end
%end
%%--------------------------------
%if model.use.segments.transexon
%  equiv = [equiv;
%           model.segments.utr5exon, model.segments.transexon] ;     
%  different = [different;
%               model.segments.intergenic, model.segments.transexon, 1] ;
%  if isfield(model.segments, 'intergenictrans'),
%    different = [different;
%                 model.segments.intergenictrans, model.segments.transexon, 1];    
%  end ;
%end
%--------------------------------

model.loss.segments = ones(model.cnt_segments) - eye(model.cnt_segments);
for i=1:size(different,1)
  model.loss.segments(different(i,1)+1, different(i,2)+1) = different(i,3) ;
  % model.loss.segments(different(i,1)+1, different(i,2)+1) = different(i,3)+1 ;warning('make loss matrix asymetric to enhance sensitivity')
  model.loss.segments(different(i,2)+1, different(i,1)+1) = different(i,3) ;
end ;

for i=1:size(equiv,1)
  model.loss.segments(equiv(i,1)+1, equiv(i,2)+1) = 0 ;
  model.loss.segments(equiv(i,2)+1, equiv(i,1)+1) = 0 ;
end ;


% tripple loss for coding
%model.loss.segments(model.segments.cds_exon+1, :) = 3*model.loss.segments(model.segments.cds_exon+1, :) ;
%model.loss.segments(:, model.segments.cds_exon+1) = 3*model.loss.segments(:, model.segments.cds_exon+1) ;


% define nucleotide loss only for utrs
if 1
	utr_nuc_loss = 1e-3;
	nuc_loss = zeros(model.cnt_segments);
	nuc_loss(model.segments.intergenic+1, model.segments.utr5exon+1) = utr_nuc_loss;
	nuc_loss(model.segments.utr5exon+1, model.segments.intergenic+1) = utr_nuc_loss;
	nuc_loss(model.segments.intergenic+1, model.segments.utr3exon+1) = utr_nuc_loss;
	nuc_loss(model.segments.utr3exon+1, model.segments.intergenic+1) = utr_nuc_loss;
	model.loss.segments = [model.loss.segments nuc_loss] ;
else
	% first part is for segments, the second for nucleotides
	if isfield(model.use, 'nucleotide_loss_factor')
	  model.loss.segments = [model.loss.segments model.use.nucleotide_loss_factor*model.loss.segments];
	else
	  model.loss.segments = [model.loss.segments 0.0*model.loss.segments] ;
	end
end

if model.use.signals.tss
  model.loss.zero_range.tss = 10 ;
end
model.loss.zero_range.tis = 0;
model.loss.zero_range.don = 0;
model.loss.zero_range.acc = 0;
if model.use.signals.trans
  model.loss.zero_range.trans = 0;
end
model.loss.zero_range.cdsStop = 0;
if model.use.signals.polya
  model.loss.zero_range.polya = 10 ;
end
if model.use.signals.cleave
  model.loss.zero_range.cleave = 10 ;
end

model.loss.confirmed = 1 ;
model.loss.confirmed_tol = 1000 ; % region around confirmed gene
model.loss.unconfirmed = 1 ;

% model.loss.mustbe_corr.intergenic = 0 ;
% model.loss.mustbe_corr.utr5exon   = 0;
% model.loss.mustbe_corr.transexon  = 0 ;
% model.loss.mustbe_corr.cds_exon   = 1 ;
% model.loss.mustbe_corr.utr3exon   = 0 ;
% model.loss.mustbe_corr.polya_tail = 0 ;
% model.loss.mustbe_corr.intercistronic  = 0 ;
% model.loss.mustbe_corr.intron     = 1 ;
return ;
