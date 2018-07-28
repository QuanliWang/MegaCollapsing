#master prorgam
library(MegaCollapsing)
sample_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/PEDMAP_v3_GGE_incCovar.txt"
matrix_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/4_v3_GGE_CCDS_syn_0005_ExACEVS0_matrix.txt"
input.data <- read.collapsing.data(sample_file,NULL,matrix_file, "IID")

meta_gene_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/CAKUT_msigdb.v6.0.symbols"
gene.sets <- read.gene.sets(meta_gene_file)

#one can simulate a random gene sets for teting
#gene.sets <-random.gene.sets(gene.sets, rownames(input.data$non.syn))

#create mega matrix (syn can be missing, but not non.syn, one can fake syn to be non.syn though)
#mega.matrix <- as.mega.matrix(gene.sets[1:20], input.data)
mega.matrix <- as.mega.matrix(gene.sets, input.data)

## save output
write.table(mega.matrix$mega.non.syn,file="collapsing_matrix.tsv",sep ="\t")
write.table(input.data$sample.list,file = "samplelist.tsv",sep = "\t",row.names = FALSE,col.names = FALSE)

mega <- read.data("samplelist.tsv", "collapsing_matrix.tsv")
ps <- get.pvalues(mega$data, mega$is.case, n.permutations = 1000)

library(QQperm)
estlambda2(ps$observed, ps$perm, plot = TRUE, filter = TRUE, adjust.xy = FALSE)


