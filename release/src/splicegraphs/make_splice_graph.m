function genes = make_splice_graph(genes)

  addpath ~/svn/projects/splicing/splicegraphs/;
  addpath ~/svn/projects/splicing/splicegraphs/detect_altsplice;
  
  build_splice_graph ;
  infer_splice_graph;
  genes= alt_const(genes) ;