#master prorgam
library(MegaCollapsing)

sample_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/PEDMAP_v3_GGE_incCovar.txt"
syn_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/4_v3_GGE_CCDS_syn_0005_ExACEVS0_matrix.txt"
non_syn_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/1_v3_GGE_CCDS_hotzone_0005_ExACEVS0_matrix.txt"
meta_gene_file <- "/Users/qw2192/Desktop/AZ/MegaCollapsing/sup/CAKUT_msigdb.v6.0.symbols"

input.data <- read.collapsing.data(sample_file, syn_file, non_syn_file, sample.column = "IID")
gene.sets <- read.gene.sets(meta_gene_file)
rm(sample_file,syn_file,non_syn_file,meta_gene_file)

## generate meta gene matrices
mega.matrix <- as.mega.matrix(gene.sets, input.data)


## save output
write.table(res,file="syn_collapsing_matrix2.tsv",sep ="\t")
write.table(sample_list,file = "samplelist.tsv",sep = "\t",row.names = FALSE,col.names = FALSE)
write.table(res,file="non_syn_collapsing_matrix2.tsv",sep ="\t")

save.image(file = "all.data.RData")


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
