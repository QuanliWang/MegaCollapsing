#master prorgam

library(MetaCollapsing)

sample_file <- "/Users/qw2192/Desktop/MetaCollapsing/PEDMAP_v3_GGE_incCovar.txt"
sample_list <- read.sample.list(sample_file)

#matrix_file <- "/Users/qw2192/Desktop/MetaCollapsing/1_v3_GGE_CCDS_hotzone_0005_ExACEVS0_matrix.txt"
matrix_file <- "/Users/qw2192/Desktop/MetaCollapsing/4_v3_GGE_CCDS_syn_0005_ExACEVS0_matrix.txt"
collapsing_matrix <- read.collapsing.matrix(matrix_file)

#reorder columns
common_sample_list <-intersect(colnames(collapsing_matrix),sample_list[,"IID"])
collapsing_matrix <- collapsing_matrix[,common_sample_list]

#redorder rows
rownames(sample_list) <- sample_list[,"IID"]
sample_list <- sample_list[common_sample_list,]

meta_gene_file <- "/Users/qw2192/Desktop/MetaCollapsing/CAKUT_msigdb.v6.0.symbols"
gene_sets <- read.gene.sets(meta_gene_file)

#gene_sets_backup <- gene_sets

## random sample some gene set
#gene_set_sizes <- unlist(lapply(gene_sets,length))
#all_genes <- rownames(collapsing_matrix)
#for (i in 1:length(gene_sets)) {
#  gene_sets[[i]][3:gene_set_sizes[i]] <- sample(all_genes,gene_set_sizes[i]-2)
#}

### code test
#all_genes <- rownames(collapsing_matrix)
#for (i in 1:length(gene_sets)) {
#  gene_sets[[i]] <- gene_sets[[i]][1:3]
#  gene_sets[[i]][3] <- sample(all_genes,1)
#}

## renerate meta gene matrices
Meta_matrix <- matrix(NA,nrow = length(gene_sets), ncol = length(common_sample_list))
Meta_matrix2 <- matrix(NA,nrow = length(gene_sets), ncol = length(common_sample_list))
meta_gene_names<- sapply(gene_sets, function(x) x[1])

missing_genes = list()
for (i in 1:length(gene_sets)) {
  print(i)
  genes <- gene_sets[[i]][3:length(gene_sets[[i]])]
  common_genes <- intersect(genes,rownames(collapsing_matrix))
  missing_genes[[meta_gene_names[i]]] <- setdiff(genes,rownames(collapsing_matrix))
  if (length(common_genes)>0) {
    temp <- collapsing_matrix[common_genes,]
    Meta_matrix[i,] <- colSums(temp)
    temp[temp > 0] <- 1
    Meta_matrix2[i,] <- colSums(temp)
  }
}

## save output
res <- as.data.frame(Meta_matrix)
colnames(res) <- common_sample_list
rownames(res) <- meta_gene_names
write.table(res,file="collapsing_matrix1.tsv",sep ="\t")

res <- as.data.frame(Meta_matrix2)
colnames(res) <- common_sample_list
rownames(res) <- meta_gene_names
write.table(res,file="collapsing_matrix2.tsv",sep ="\t")

write.table(sample_list,file = "samplelist.tsv",sep = "\t",row.names = FALSE,col.names = FALSE)

igm.data <- read.data("samplelist1.tsv", "collapsing_matrix1.tsv")
na.row <- rowSums(is.na(igm.data$data))

to_be_removed <- which(na.row == dim(igm.data$data)[2])
if (length(to_be_removed) > 0) {
  igm.data$data <- igm.data$data[-to_be_removed,]
}

ps <- get.pvalues(igm.data$data, igm.data$is.case, n.permutations = 1000)
#save.image(file = "last1.RData")


library(QQperm)
estlambda2(ps$observed, ps$perm, plot = TRUE, filter = TRUE, adjust.xy = FALSE)


