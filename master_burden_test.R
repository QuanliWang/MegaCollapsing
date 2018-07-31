#master prorgam
library(MegaCollapsing)

#sample_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/PEDMAP_v3_GGE_incCovar.txt"
#syn_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/4_v3_GGE_CCDS_syn_0005_ExACEVS0_matrix.txt"
#non_syn_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/1_v3_GGE_CCDS_hotzone_0005_ExACEVS0_matrix.txt"
#mega_gene_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/CAKUT_msigdb.v6.0.symbols"
#exclude_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/ExcludeGenes.txt"

sample_file <- "/Users/qw2192/Desktop/new_study/PEDMAP_SAGA_incCovar.txt"
syn_file <- "/Users/qw2192/Desktop/new_study/2018-07-11_09-41-31_dominantSynonymous_matrix.txt"
non_syn_file <- "/Users/qw2192/Desktop/new_study/2018-07-11_08-39-46_dominantFlexiblePolyphenDamaging_matrix.txt"
mega_gene_file <- "/Users/qw2192/Desktop/new_study/pruned_msigdb.v6.0.symbols"
exclude_file <- NULL

input.data <- read.collapsing.data(sample_file, syn_file, non_syn_file, sample.column = "IID")
if(!is.null(exclude_file)) {
  input.data <- exclude.genes(input.data, exclude_file)
}
gene.sets <- read.gene.sets(mega_gene_file)
rm(sample_file,syn_file,non_syn_file,mega_gene_file,exclude_file)

## generate meta gene matrices
mega.matrix <- as.mega.matrix(gene.sets, input.data)

## save output
#write.table(input.data$sample.list, file = "samplelist.tsv",sep = "\t",row.names = FALSE,col.names = FALSE)
#write.table(mega.matrix$mega.syn,file="syn_collapsing_matrix.tsv",sep ="\t")
#write.table(mega.matrix$mega.non.syn ,file="non_syn_collapsing_matrix.tsv",sep ="\t")

#just save a copy of the data, as it took so long to compute
#save.image(file = "all.data.RData")

# fit regression models
X <- c("Gender")
Y = NULL #default to sample_list[,'Status'] - 1
result <- burden.test(input.data$sample.list, gene.sets[1:10], mega.matrix$mega.syn,mega.matrix$mega.non.syn, X, Y)

write.table(result$p_values,file="p_values_from_mega_collapsing.tsv",sep ="\t")



