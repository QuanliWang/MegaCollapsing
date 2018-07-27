#master prorgam

library(MegaCollapsing)

sample_file <- "/Users/qw2192/Desktop/MetaCollapsing/sup/PEDMAP_v3_GGE_incCovar.txt"
sample_list <- read.sample.list(sample_file)

matrix_file <- "/Users/qw2192/Desktop/MetaCollapsing/sup/1_v3_GGE_CCDS_hotzone_0005_ExACEVS0_matrix.txt"
non_syn_collapsing_matrix <- read.collapsing.matrix(matrix_file)

matrix_file <- "/Users/qw2192/Desktop/MetaCollapsing/sup/4_v3_GGE_CCDS_syn_0005_ExACEVS0_matrix.txt"
syn_collapsing_matrix <- read.collapsing.matrix(matrix_file)

#reorder columns
common_sample_list <-intersect(colnames(non_syn_collapsing_matrix),sample_list[,"IID"])
non_syn_collapsing_matrix <- non_syn_collapsing_matrix[,common_sample_list]
syn_collapsing_matrix <- syn_collapsing_matrix[,common_sample_list]

#redorder rows
rownames(sample_list) <- sample_list[,"IID"]
sample_list <- sample_list[common_sample_list,]

meta_gene_file <- "/Users/qw2192/Desktop/MetaCollapsing/sup/CAKUT_msigdb.v6.0.symbols"
gene_sets <- read.gene.sets(meta_gene_file)


## renerate meta gene matrices
syn_Meta_matrix <- matrix(NA,nrow = length(gene_sets), ncol = length(common_sample_list))
syn_Meta_matrix2 <- matrix(NA,nrow = length(gene_sets), ncol = length(common_sample_list))
meta_gene_names<- sapply(gene_sets, function(x) x[1])

non_syn_Meta_matrix <- matrix(NA,nrow = length(gene_sets), ncol = length(common_sample_list))
non_syn_Meta_matrix2 <- matrix(NA,nrow = length(gene_sets), ncol = length(common_sample_list))

missing_genes = list()
for (i in 1:length(gene_sets)) {
  print(i)
  genes <- gene_sets[[i]][3:length(gene_sets[[i]])]
  common_genes <- intersect(genes,rownames(syn_collapsing_matrix))
  missing_genes[[meta_gene_names[i]]] <- setdiff(genes,rownames(syn_collapsing_matrix))
  if (length(common_genes)>0) {
    temp <- syn_collapsing_matrix[common_genes,]
    syn_Meta_matrix[i,] <- colSums(temp)
    temp[temp > 0] <- 1
    syn_Meta_matrix2[i,] <- colSums(temp)

    temp <- non_syn_collapsing_matrix[common_genes,]
    non_syn_Meta_matrix[i,] <- colSums(temp)
    temp[temp > 0] <- 1
    non_syn_Meta_matrix2[i,] <- colSums(temp)
  }
}


## save output
res <- as.data.frame(syn_Meta_matrix)
colnames(res) <- common_sample_list
rownames(res) <- meta_gene_names
write.table(res,file="syn_collapsing_matrix1.tsv",sep ="\t")

res <- as.data.frame(syn_Meta_matrix2)
colnames(res) <- common_sample_list
rownames(res) <- meta_gene_names
write.table(res,file="syn_collapsing_matrix2.tsv",sep ="\t")

write.table(sample_list,file = "samplelist.tsv",sep = "\t",row.names = FALSE,col.names = FALSE)

## save output
res <- as.data.frame(non_syn_Meta_matrix)
colnames(res) <- common_sample_list
rownames(res) <- meta_gene_names
write.table(res,file="non_syn_collapsing_matrix1.tsv",sep ="\t")

res <- as.data.frame(non_syn_Meta_matrix2)
colnames(res) <- common_sample_list
rownames(res) <- meta_gene_names
write.table(res,file="non_syn_collapsing_matrix2.tsv",sep ="\t")


# use matrix2 for further analysis, as this is the counts for number of genes with variants

rm(non_syn_Meta_matrix,syn_Meta_matrix,missing_genes,res,temp)

# fit regression models
Y <- sample_list[,'Status'] - 1
p_values = matrix(NA,length(gene_sets),4)
rownames(p_values) <- meta_gene_names
colnames(p_values) <- c('(Intercept)','V1','gender','sV1') #a hack
gender = sample_list[,'Gender']
for (i in 1:length(gene_sets)) {
  print(i)
  V1 = non_syn_Meta_matrix2[i,]
  sV1 = syn_Meta_matrix2[i,]
  try({
  V1_glm <- glm( Y ~ V1 + gender + sV1, family = binomial(link = logit))
  coef <- summary(V1_glm)$coef
  p_values[i,rownames(coef)] <- coef[,4]
  }, silent=TRUE)
}

write.table(p_values,file="p_values_from_collapsing_matrix2.tsv",sep ="\t")
