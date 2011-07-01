function gene_eval(anno_dir, pred_dir)

l = load([anno_dir '/genes.mat'], 'genes');
ll = load([pred_dir '/genes.mat'], 'genes');

[correct SN SP F] = eval_fast(ll.genes, l.genes);
