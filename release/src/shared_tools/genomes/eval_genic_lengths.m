function lengths = eval_genic_lengths(genes)
  tmp=length(genes);
  lengths.utr5exon=zeros(1,tmp*2);
  utr5it=1;
  lengths.single_cds_exon=zeros(1,tmp);
  singleit=1;
  lengths.first_cds_exon=zeros(1,tmp);
  firstit=1;
  lengths.middle_cds_exon=zeros(1,tmp*10);
  middelit=1;
  lengths.last_cds_exon=zeros(1,tmp);
  lastit=1;
  lengths.utr3exon=zeros(1,tmp*2);
  utr3it=1;
  lengths.intron=zeros(1,tmp*10)
  intronit=1;
  for i=1:length(genes)
    for j=1:length(genes(i).transcripts)
      for k=1:getcolumn(size(genes(i).utr5_exons{j}),1)
	lengths.utr5exon(utr5it)=genes(i).utr5_exons{j}(k,2)-genes(i).utr5_exons{j}(k,1);
        utr5it=utr5it+1;
      end
      num_cds_exons = getcolumn(size(genes(i).cds_exons{j}),1);
      if num_cds_exons==1
	lengths.single_cds_exon(singleit)=genes(i).cds_exons{j}(1,2)-genes(i).cds_exons{j}(1,1);
	singleit=singleit+1;
      else
	lengths.first_cds_exon(firstit)=genes(i).cds_exons{j}(1,2)-genes(i).cds_exons{j}(1,1);
	firstit=firstit+1;
	for k=2:num_cds_exons-1
	  lengths.middle_cds_exon(middleit)=genes(i).cds__exons{j}(k,2)-genes(i).cds__exons{j}(k,1);
          utr5it=utr5it+1;
	end
	lengths.last_cds_exon(firstit)=genes(i).cds_exons{j}(end,2)-genes(i).cds_exons{j}(end,1);
        lastit=lastit+1;
      end
      for k=1:getcolumn(size(genes(i).utr3_exons{j}),1)
	lengths.utr3exons(utr3it)=genes(i).utr3_exons{j}(k,2)-genes(i).utr3_exons{j}(k,1);
        utr3it=utr3it+1;
      end
      for k=1:getcolumn(size(genes(i).exons{j}),1)-1
	lengths.intron(intronit)=genes(i).exons{j}(k+1,1)-genes(i).utr3_exons{j}(k,2);
        intronit=intronit+1;
      end
    end
  end
  % cut away unused space
  lengths.utr5exon=lengths.utr5exon(1:utr5it-1);
  lengths.single_cds_exon=lengths.single_cds_exon(1:singleit-1);
  lengths.first_cds_exon=lengths.first_cds_exon(1:firstit-1);
  lengths.middle_cds_exon=lengths.middle_cds_exon(1:middelit-1);
  lengths.last_cds_exon=lengths.last_cds_exon(1:lastit-1);
  lengths.utr3exon=lengths.utr3exon(1:utr3it-1);
  lengths.intron=lengths.intron(1:intronit-1);
end