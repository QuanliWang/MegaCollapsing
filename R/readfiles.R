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
  if (length(offending_sample_name_index) > 0) {
    samplelist[,"IID"][offending_sample_name_index] <-
      sapply(samplelist[offending_sample_name_index,"IID"], function(x) paste("X",x,sep=""))
  }
  return(samplelist)
}

read.collapsing.matrix <- function(matrix_file) {
  mat <- read.table(matrix_file,header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  rownames(mat) <- mat[,"sample.gene"]
  return(mat)
}

read.collapsing.data <- function(samples, syn, non.syn, sample.column = "IID") {
  if (file.exists(samples)) {
    sample_list <- read.sample.list(samples)
  } else {
    stop("input sample file not found")
  }

  if (!is.element(sample.column, colnames(sample_list))) {
    stop("sample ID column not found")
  }

  if (!is.null(syn)) {
    if (file.exists(syn)) {
      syn.mat <- read.collapsing.matrix(syn)
    } else {
      stop("input syn matrix file not found")
    }
  }

  if (file.exists(non.syn)) {
    non.syn.mat <- read.collapsing.matrix(non.syn)
  } else {
    stop("input non.syn matrix file not found")
  }

  if (!is.null(syn)) {
    if (!all.equal(dim(syn.mat),dim(non.syn.mat))) {
      stop("syn and non.syn do not have matching dimensions")
    }
  }

  #reorder columns
  common_sample_list <-intersect(colnames(non.syn.mat),sample_list[,sample.column])
  if (!is.null(syn)) {
    syn.mat <- syn.mat[,common_sample_list]
  }
  non.syn.mat <- non.syn.mat[,common_sample_list]

  #redorder rows
  rownames(sample_list) <- sample_list[,sample.column]
  sample_list <- sample_list[common_sample_list,]
  if (is.null(syn)) {
    return(list( sample.list = sample_list, non.syn = non.syn.mat))
  } else {
    return(list( sample.list = sample_list, syn = syn.mat, non.syn = non.syn.mat))
  }
}
