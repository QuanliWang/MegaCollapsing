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

  na.row <- rowSums(is.na(matrix))
  to_be_removed <- which(na.row == dim(matrix)[2])
  if (length(to_be_removed) > 0) {
    matrix <- matrix[-to_be_removed,]
  }

  out <- list()
  out$data <- matrix
  out$is.case <- is.case


  out
}

get.pvalues <-function(matrix, is.case, n.permutations = 1000, with_contingency = FALSE) {
  #Number of cases and controls
  matrix[matrix>0] <- 1
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

  total.1 <- rowSums(matrix)
  if (n.permutations > 0) {
  #permutation, save all p-values just in case median will be needed later on
    P.Values <- matrix(1,dim(matrix)[1],n.permutations)

    for (i in 1: n.permutations) {
      K <- sample.int(n.samples, size = n.cases, replace = FALSE)
      Labels.1.1 <- rowSums(matrix[,K])
      Labels.0.1 <- total.1 - Labels.1.1
      P.Values[,i] <- sort(Fisher.precompute[cbind(Labels.1.1+1,Labels.0.1+1)])
    }
    P.perm <- rowMeans(P.Values)
    #P.perm <- sapply(P.Values, function(x) quantile(x,c(0.025, 0.50, 0.975)))
  }


  #compute observed (true case control configration) p-values
  K <- which(is.case)
  Labels.1.1 <- rowSums(matrix[,K])
  Labels.0.1 <- total.1 - Labels.1.1

  if (n.permutations > 0) {
    P.observed <- sort(Fisher.precompute[cbind(Labels.1.1+1,Labels.0.1+1)])
  } else {
    P.observed <- Fisher.precompute[cbind(Labels.1.1+1,Labels.0.1+1)] #no sorting
  }
  names(P.observed) <- rownames(matrix)

  if (with_contingency) {
    contingency_table <- matrix(NA, nrow = length(P.observed),4)
    contingency_table[,1] <- Labels.1.1
    contingency_table[,2] <- n.cases - contingency_table[,1]
    contingency_table[,3] <- Labels.0.1
    contingency_table[,4] <- n.controls - contingency_table[,3]
    rownames(contingency_table) <- rownames(matrix)
  }

  out <- list()
  if (n.permutations > 0) {
    out$perm <- P.perm
    #out$perm <- P.perm[,2]
    #out$perm.LCI <- P.perm[,1]
    #out$perm.UCI <- P.perm[,3]
  }
  out$observed <- P.observed
  if (with_contingency) {
    out$contingency.table <- contingency_table
  }
  out
}

## random sample some gene set
random.gene.sets <- function(gene.sets, genes) {
  gene.set.sizes <- unlist(lapply(gene.sets,length))
  for (i in 1:length(gene.sets)) {
    gene.sets[[i]][3:gene.set.sizes[i]] <- sample(genes,gene.set.sizes[i]-2)
  }
  return(gene.sets)
}
