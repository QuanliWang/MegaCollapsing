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
