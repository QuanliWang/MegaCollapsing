#master prorgam

#library(MetaCollapsing)
read.seed.list <- function(seed_set) {
  # tab delimited but read it as csv file
  seedgenes <-  read.csv(seed_set,header = FALSE, quote = "",stringsAsFactors = FALSE)
  seedgenes <- lapply(seedgenes,as.character)[[1]]
  ret = list()
  for (i in 1:length(seedgenes)) {
    ret[[i]] <- unlist(strsplit(gsub(" ", "",seedgenes[i]),"/"))
  }
  return(ret)
}

read.gene.sets <- function(meta_gene_file) {
  # tab delimited but read it as csv file
  genesets <-  read.csv(meta_gene_file,header = FALSE, stringsAsFactors = FALSE)
  genesets <- genesets[,1]
  genesets <- strsplit(genesets,'\t',fixed = TRUE)
  return(genesets)
  #save(genesets, file = "genesets.RData")
}

read.sample.list <- function(sample_file) {
  samplelist <-  read.table(sample_file,header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  offending_sample_name_index <- grep("^[0-9]",samplelist[,"IID"])
  samplelist[,"IID"][offending_sample_name_index] <-
    sapply(samplelist[offending_sample_name_index,"IID"], function(x) paste("X",x,sep=""))
  return(samplelist)
}

read.collapsing.matrix <- function(matrix_file) {
  mat <- read.table(matrix_file,header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  rownames(mat) <- mat[,"sample.gene"]
  return(mat)
}

read.data <- function(sample.file, matrix.file, filter.list = NULL) {
  #get sample and phenotype label information
  ped <- read.table(sample.file,header = FALSE, sep = '\t', as.is = TRUE)
  samples <- ped[,2]
  is.case <- as.matrix(ped[,6] ==2)
  rownames(is.case) <- samples

  #get count matrix
  data <- read.table(matrix.file,header = TRUE, sep = '\t', as.is = TRUE)
  genes <- data[,1]
  data <- data[,samples]
  data[data > 0] <- 1

  #convert data to matrix
  matrix <- as.matrix(data[])
  rownames(matrix)<-genes

  if (!is.null(filter.list)) {
    genes.common <- intersect(genes,filter.list)
    matrix <- matrix[genes.common,]
  }
  out <- list()
  out$data <- matrix
  out$is.case <- is.case
  out
}

get.pvalues <-function(matrix, is.case, n.permutations = 1000) {
  #Number of cases and controls
  n.samples <- length(is.case)
  n.cases <- sum(is.case)
  n.controls <- n.samples - n.cases

  #pre-compute all possible contingency.table for performance
  #this will create a look up table
  contingency.table = matrix(0,2,2)
  max.counts <- max(rowSums(matrix))
  Fisher.precompute <- matrix(0, min(n.cases, max.counts) + 1,
                              min(n.controls, max.counts) + 1)
  for (i in 1: dim(Fisher.precompute)[1]) {
    for (j in 1:dim(Fisher.precompute)[2]) {
      contingency.table[1,1] <- i-1
      contingency.table[1,2] <- n.cases - contingency.table[1,1]
      contingency.table[2,1] <- j -1
      contingency.table[2,2] <- n.controls - contingency.table[2,1]
      #if (all(contingency.table >=0)) {
      Fisher.precompute[i,j] <- fisher.test(contingency.table)$p.value
      #}
    }
  }

  #permutation, save all p-values just in case median will be needed later on
  P.Values <- matrix(1,dim(matrix)[1],n.permutations)
  total.1 <- rowSums(matrix)
  for (i in 1: n.permutations) {
    K <- sample.int(n.samples, size = n.cases, replace = FALSE)
    Labels.1.1 <- rowSums(matrix[,K])
    Labels.0.1 <- total.1 - Labels.1.1
    P.Values[,i] <- sort(Fisher.precompute[cbind(Labels.1.1+1,Labels.0.1+1)])
  }
  P.perm <- rowMeans(P.Values)
  #P.perm <- sapply(P.Values, function(x) quantile(x,c(0.025, 0.50, 0.975)))


  #compute observed (true case control configration) p-values
  K <- which(is.case)
  Labels.1.1 <- rowSums(matrix[,K])
  Labels.0.1 <- total.1 - Labels.1.1
  P.observed <- sort(Fisher.precompute[cbind(Labels.1.1+1,Labels.0.1+1)])

  out <- list()
  out$perm <- P.perm
  #out$perm <- P.perm[,2]
  #out$perm.LCI <- P.perm[,1]
  #out$perm.UCI <- P.perm[,3]
  out$observed <- P.observed
  out
}


sample_file <- "/Users/qw2192/Desktop/MetaCollapsing/sup/PEDMAP_v3_GGE_incCovar.txt"
sample_list <- read.sample.list(sample_file)

#matrix_file <- "/Users/qw2192/Desktop/MetaCollapsing/sup/1_v3_GGE_CCDS_hotzone_0005_ExACEVS0_matrix.txt"
matrix_file <- "/Users/qw2192/Desktop/MetaCollapsing/sup/4_v3_GGE_CCDS_syn_0005_ExACEVS0_matrix.txt"
collapsing_matrix <- read.collapsing.matrix(matrix_file)

#reorder columns
common_sample_list <-intersect(colnames(collapsing_matrix),sample_list[,"IID"])
collapsing_matrix <- collapsing_matrix[,common_sample_list]

#redorder rows
rownames(sample_list) <- sample_list[,"IID"]
sample_list <- sample_list[common_sample_list,]

meta_gene_file <- "/Users/qw2192/Desktop/MetaCollapsing/sup/CAKUT_msigdb.v6.0.symbols"
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

igm.data <- read.data("samplelist.tsv", "collapsing_matrix1.tsv")
na.row <- rowSums(is.na(igm.data$data))

to_be_removed <- which(na.row == dim(igm.data$data)[2])
if (length(to_be_removed) > 0) {
  igm.data$data <- igm.data$data[-to_be_removed,]
}

ps <- get.pvalues(igm.data$data, igm.data$is.case, n.permutations = 1000)
#save.image(file = "last1.RData")


library(QQperm)
estlambda2(ps$observed, ps$perm, plot = TRUE, filter = TRUE, adjust.xy = FALSE)


