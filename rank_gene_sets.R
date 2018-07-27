#master prorgam
meta_gene_file <- "/Users/qw2192/Desktop/MetaCollapsing/CAKUT_msigdb.v6.0.symbols"
gene_sets <- read.gene.sets(meta_gene_file)
meta_gene_names<- sapply(gene_sets, function(x) x[1])

seed_set <- "/Users/qw2192/Desktop/MetaCollapsing/trial_seed_n150.txt"
seed_genes <- read.seed.list(seed_set)


#split seed genes into two groups
syn_count <- sapply(seed_genes, length)
seed_set_one <- unlist(seed_genes[syn_count == 1])
seed_set_two <- seed_genes[syn_count > 1]

seed_set_flat <- unlist(seed_genes)

out_matrix <- matrix(NA,nrow = length(gene_sets), ncol = 2)

#compute measure 1: Jaccard Index
for (i in 1:length(gene_sets)) {
  print(i)
  genes <- gene_sets[[i]][3:length(gene_sets[[i]])]
  si = length(intersect(genes,seed_set_flat))

  #a bit more complicated to compute unions
  su <- length(union(genes,seed_set_one))
  for (j in 1:length(seed_set_two)) {
    if (length(intersect(genes,seed_set_two[[j]]))==0) {
      su <- su + 1
    }
  }
  out_matrix[i,1] <- si / su
}

#compute measure 2: Petrovski Index
# compute the size of seed_genes
sb <- length(seed_genes)
for (i in 1:length(gene_sets)) {
  print(i)
  genes <- gene_sets[[i]][3:length(gene_sets[[i]])]
  sa <- length(genes)
  saminusb <- sa - length(intersect(genes,seed_set_flat))

  #a bit more complicated to compute sbminusa
  sbminusa <- sb - length(intersect(seed_set_one,genes))
  for (j in 1:length(seed_set_two)) {
    if (length(intersect(genes,seed_set_two[[j]]))>0) {
      sbminusa <- sbminusa - 1
    }
  }
  out_matrix[i,2] <-  1-(saminusb / sa + sbminusa/sb) / 2
}

## save output
res <- as.data.frame(out_matrix)
rownames(res) <- meta_gene_names
res <- cbind(rownames(res),res)
colnames(res) <- c("gene_set","Jaccard_Index", "Petrovski_Index")
write.table(res,file="gene_set_ranks.csv",sep =",",row.names = FALSE)





