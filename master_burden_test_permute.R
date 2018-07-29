library(MegaCollapsing)

sample_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/PEDMAP_v3_GGE_incCovar.txt"
syn_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/4_v3_GGE_CCDS_syn_0005_ExACEVS0_matrix.txt"
non_syn_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/1_v3_GGE_CCDS_hotzone_0005_ExACEVS0_matrix.txt"
mega_gene_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/CAKUT_msigdb.v6.0.symbols"
exclude_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/ExcludeGenes.txt"
#exclude_file <- NULL

input.data <- read.collapsing.data(sample_file, syn_file, non_syn_file, sample.column = "IID")
if(!is.null(exclude_file)) {
  input.data <- exclude.genes(input.data, exclude_file)
}
gene.sets <- read.gene.sets(mega_gene_file)
rm(sample_file,syn_file,non_syn_file,mega_gene_file,exclude_file)

## generate meta gene matrices
mega.matrix <- as.mega.matrix(gene.sets, input.data)

# fit regression models
X <- c("Gender")
Y = NULL #default to sample_list[,'Status'] - 1

# fit regression models
perm.result <- burden.test.with.permutation(input.data$sample.list, gene.sets[1:10],
                                            mega.matrix$mega.syn[1:10,], mega.matrix$mega.non.syn[1:10,],
                                            X, Y, n.permutations = 100)
write.table(perm.result$p_summary,file="p_perm_null_100.tsv",sep ="\t")

write.table(perm.result$matrix,file="p_mat.tsv",sep ="\t")

#now get a list gene sets that has p-values less than 0.001
Index <- which(perm.result$observed <= 0.001)
if (length(Index) > 0) {
  gs001 <- gene.sets[Index]
  non_syn_trimmed <- mega.matrix$mega.non.syn[Index,]
  syn_trimmed <- mega.matrix$mega.syn[Index,]
  n.permutations <- 1000

  perm.result.gs001 <- burden.test.with.permutation(input.data$sample.list, gs001,
                                                    syn_trimmed, non_syn_trimmed,
                                                    X, Y, n.permutations = n.permutations)

  write.table(perm.result$matrix,file="p_perm_0.001.tsv",sep ="\t")
}
