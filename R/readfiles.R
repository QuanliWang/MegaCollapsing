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

.read.gene.sets2 <- function(mega_gene_file) {
  # tab delimited but read it as csv file
  genesets <-  read.csv(mega_gene_file,header = FALSE, stringsAsFactors = FALSE)
  genesets <- genesets[,1]
  genesets <- strsplit(genesets,'\t',fixed = TRUE)
  return(genesets)
}

read.gene.sets <- function(mega_gene_file) {
  geneset.list <- .read.gene.sets2(mega_gene_file)
  gene.set.names <- lapply(geneset.list, function(x) x[1])
  gene.sets <- lapply(geneset.list, function(x) x[3:length(x)])
  names(gene.sets) <- gene.set.names
  unique.genes <- sort(unique(unlist(gene.sets)))
  a <- lapply(gene.sets, function(x,y) {match(x,y)}, y = unique.genes)

  f <- function(x, n) {
      r <- rep(FALSE,1,n)
      r[x] <- TRUE
      r
  }
  mat <- t(sapply(a, f, n = length(unique.genes)))
  mat <- as(mat, "sparseMatrix")
  dimnames(mat) <- list(gene.set.names,unique.genes)

  return(list(gene.sets = gene.sets, mat = mat))
}

read.sample.list <- function(sample_file, sample_ID, fix.sample.names = TRUE) {
  samplelist <-  read.table(sample_file,header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  if (!sample_ID %in% names(samplelist) && length(names(samplelist)) == 8) {
    samplelist <-  read.table(sample_file,header = FALSE, stringsAsFactors = FALSE, sep = "\t")
    names(samplelist) <- c("FID",sample_ID,"PID", "MID", "sex", "status","Seq","Kit")
  }
  if (fix.sample.names) {
      offending_sample_name_index <- grep("^[0-9]",samplelist[,sample_ID])
      samplelist[,sample_ID][offending_sample_name_index] <-
        sapply(samplelist[offending_sample_name_index,sample_ID], function(x) paste("X",x,sep=""))
  }
  return(samplelist)
}

read.collapsing.matrix <- function(matrix_file, trunk.size = 10000) {
  ## figure out the number of columns/samples
  mat <- fread(matrix_file, na.strings = c("NA",""), stringsAsFactors = FALSE, sep = "\t", header = TRUE, nrows = 1)
  samples <- names(mat)[-1]
  number.columns <- dim(mat)[2]

  #read in gene names
  genes <- fread(matrix_file, na.strings = c("NA",""), stringsAsFactors = FALSE, sep = "\t", header = TRUE, select = 1)
  genes <- unlist(genes)

  result <- NULL
  number.of.trunks <- (number.columns - 1) %/% trunk.size
  if (number.of.trunks * trunk.size < number.columns - 1) {
      number.of.trunks <- number.of.trunks + 1
  }
  for (i in 1:number.of.trunks) {
      trunk.start <- (i - 1) * trunk.size + 2
      trunk.end <- i * trunk.size + 1
      if (trunk.end >  number.columns) {
          trunk.end <- number.columns
      }
      if (trunk.end >= trunk.start) {
          mat <- fread(matrix_file, na.strings = c("NA",""), stringsAsFactors = FALSE, sep = "\t",
                   header = TRUE, select =  trunk.start: trunk.end)
          mat <- as(as.matrix(mat), "sparseMatrix")
          result <- cbind(result,mat)
      }
  }
  dimnames(result) <- list(genes,samples)
  return(result)
}

read.collapsing.data <- function(samples, syn, non.syn, sample.column = "IID") {
  if (file.exists(samples)) {
    sample_list <- read.sample.list(samples, sample.column, FALSE)
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
