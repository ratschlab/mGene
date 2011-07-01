function ALL = initialize_allconfs(POS,info_names)

ALL.POS = POS;  
  
ALL.GENE_ID = zeros(1,length(POS)); 
ALL.ALTGENIC = zeros(1,length(POS)); 
ALL.COV1 = zeros(length(info_names),length(POS)); 
ALL.COV2 = zeros(length(info_names),length(POS)); 
ALL.LABEL_2 =  -ones(3,length(POS)); 
ALL.ISSPLICE = zeros(1,length(POS)); 

ALL.TRANS = zeros(1,length(POS)); 
ALL.TSS = zeros(1,length(POS)); 
ALL.HAS_SIG = zeros(1,length(POS));
ALL.DIST = zeros(1,length(POS)); 
ALL.CONF_TRANS = zeros(1,length(POS)); 
