exclude.samples.by.genes <- function(input.data, exclude.file) {
  if (!file.exists(exclude.file)) {
    stop("exlude gene file not found")
  }
  exclude.list  <-  read.csv(exclude.file,header = FALSE, stringsAsFactors = FALSE)
  genes <- rownames(input.data$non.syn)
  common.genes <- intersect(genes, exclude.list[,1])
  if (length(common.genes) > 0) {
    index <- match(common.genes, genes)
    trimmed.data <- input.data$non.syn[index,]
    sample.to.exclude <- which(colSums(trimmed.data) > 0)
    if (length(sample.to.exclude) > 0) {
      input.data$sample.list <- input.data$sample.list[-sample.to.exclude,]
      input.data$syn <- input.data$syn[,-sample.to.exclude]
      input.data$non.syn <- input.data$non.syn[,-sample.to.exclude]
    }
  }
  return(input.data)
}

get.samples.by.exclude.genes <- function(non.syn, exclude) {
    genes <- rownames(non.syn)
    samples <- colnames(non.syn)
    exclude <- exclude %>% mutate(samples = "")
    for (i in 1:dim(exclude)[1]) {
        gene.list <- unlist(strsplit(exclude$genolist[i],","))
        common.genes <- intersect(genes, gene.list)
        if (length(common.genes) > 0) {
            index <- match(common.genes, genes)
            trimmed.data <- non.syn[index,]
            if (is.null(dim(trimmed.data))) {
                dim(trimmed.data) <- c(1,length(trimmed.data))
            }
            sample.to.exclude.index <- which(colSums(trimmed.data) > 0)
            if (length(sample.to.exclude.index) > 0) {
                sample.to.exclude <- samples[sample.to.exclude.index]
                exclude$samples[i] <- paste0(sample.to.exclude, collapse = "|")
            }
        }
    }
    return(exclude)
}

as.mega.matrix<- function(gene.sets, input.data) {
  sample.list <- colnames(input.data$non.syn)
  #mega.gene.names <- rownames(gene.sets$mat)

  gene.names.from.gene.sets <- colnames(gene.sets$mat)
  gene.names.from.collapsing.matrix <- rownames(input.data$non.syn)

  has.syn <- "syn" %in% names(input.data)

  common.genes <- intersect(gene.names.from.gene.sets, gene.names.from.collapsing.matrix)
  if (length(common.genes) == 0) {
      stop("no common genes were found among collapsing matrix and gene sets.")
  }

  gene.index <- match(common.genes, gene.names.from.collapsing.matrix)
  if (has.syn) {
      input.data$syn <- input.data$syn[gene.index,]
  }
  input.data$non.syn <- input.data$non.syn[gene.index,]

  gene.index <- match(common.genes, gene.names.from.gene.sets)
  gene.sets$mat <- gene.sets$mat[,gene.index]
  gene.sets$mat[gene.sets$mat] <-1
  gene.sets$mat[!gene.sets$mat] <-0
  if (has.syn) {
      mega.syn <- gene.sets$mat %*% input.data$syn
  }
  mega.non.syn <- gene.sets$mat %*% input.data$non.syn

  if (has.syn) {
    return(list(mega.syn = mega.syn, mega.non.syn = mega.non.syn, sample.list = sample.list))
  } else {
    return(list(mega.non.syn = mega.non.syn, sample.list = sample.list))
  }
}

