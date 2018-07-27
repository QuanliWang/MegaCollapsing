library(MegaCollapsing)
syn_file <- "/Users/qw2192/Desktop/AZ/mega-matrices/syn_logit.txt"
non_syn_file <- "/Users/qw2192/Desktop/AZ/mega-matrices/logit.txt"

syn <- read.table(syn_file,header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = "\t")
non_syn <- read.table(non_syn_file,header = TRUE,row.names = 1, stringsAsFactors = FALSE, sep = "\t")

syn = t(syn)
non_syn = t(non_syn)

sample_file <- "/Users/qw2192/Desktop/AZ/mega-matrices/samplelist_dom_composite_matrix.tsv"
samples <- read.table(sample_file,header = TRUE,stringsAsFactors = FALSE, sep = "\t")
rownames(samples) <- samples[,"IID"]
samples <- samples[colnames(syn),]

#setdiff(rownames(samples),colnames(syn))

meta_gene_file <- "/Users/qw2192/Desktop/AZ/Quanli/pruned_msigdb.v6.0.symbols"
gene_sets <- read.gene.sets(meta_gene_file)


# fit regression models
Y <- samples[,'Status'] - 1
p_values = matrix(NA,length(gene_sets),4)
meta_gene_names<- sapply(gene_sets, function(x) x[1])
rownames(p_values) <- meta_gene_names
colnames(p_values) <- c('(Intercept)','V1','gender','sV1') #a hack
gender = samples[,'Gender']


#perm.result <- get.pvalues_logit_null(Y[1:10], gender[1:10],non_syn[1:10,],syn[1:10,], n.permutations = 2)
perm.result <- get.pvalues_logit_null(Y, gender,non_syn,syn, n.permutations = 100)

p_values1 <- matrix(NA,length(gene_sets),2)
p_values1[,1] <- perm.result$perm
p_values1[,2] <- perm.result$observed
rownames(p_values1) <-  sapply(gene_sets, function(x) x[1])
colnames(p_values1) <- c('perm','observed')
write.table(p_values1,file="p_perm_null_100.tsv",sep ="\t")


p_mat <- perm.result$matrix
rownames(p_mat) <-  sapply(gene_sets, function(x) x[1])
write.table(p_mat,file="p_mat.tsv",sep ="\t")


library(QQperm)
estlambda2(perm.result$observed, perm.result$perm, plot = TRUE, filter = TRUE, adjust.xy = FALSE)


##########################################
for (i in 1:length(gene_sets)) {
  print(i)
  V1 = non_syn[i,]
  sV1 = syn[i,]
  try({
    V1_glm <- glm( Y ~ V1 + gender + sV1, family = binomial(link = logit))
    coef <- summary(V1_glm)$coef
    p_values[i,rownames(coef)] <- coef[,4]
  }, silent=TRUE)
}

write.table(p_values,file="p_values_from_collapsing_matrix2.tsv",sep ="\t")

#now get a list gene sets that has p-values less than 0.001
Index <- which(p_values[,"V1"] <= 0.001)
gs001 <- gene_sets[Index]
non_syn_trimmed <- non_syn[Index,]
syn_trimmed <- syn[Index,]
n.permutations <- 1000

p_values_Index = matrix(NA,length(gs001),2)
for (i in 1:length(gs001)) {
  print(i)
  V1 = non_syn_trimmed[i,]
  sV1 = syn_trimmed[i,]
  perm.result <- get.pvalues_logit(Y, V1, gender, sV1, n.permutations = n.permutations)
  p_values_Index[i,1] <- perm.result$perm
  p_values_Index[i,2] <- perm.result$observed
  #p_values[Index[i],]
}

rownames(p_values_Index) <-  sapply(gs001, function(x) x[1])
colnames(p_values_Index) <- c('perm','observed')
write.table(p_values_Index,file="p_perm_0.001.tsv",sep ="\t")
