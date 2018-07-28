as.mega.matrix<- function(gene.sets, input.data) {
  mega.gene.names<- sapply(gene.sets, function(x) x[1])
  common_sample_list <- rownames(input.data$sample.list)
  mega.syn <- matrix(NA,nrow = length(gene.sets), ncol = length(common_sample_list))
  mega.non.syn <- matrix(NA,nrow = length(gene.sets), ncol = length(common_sample_list))
  missing_genes = list()
  for (i in 1:length(gene.sets)) {
    if (i %% 10 == 0) {
      print(paste("Processing ", i, " out of ", length(gene.sets), " megagenes..."))
    }
    genes <- gene.sets[[i]][3:length(gene.sets[[i]])]
    common_genes <- intersect(genes,rownames(input.data$syn))
    missing_genes[[mega.gene.names[i]]] <- setdiff(genes,rownames(input.data$syn))
    if (length(common_genes)>0) {
      temp <- input.data$syn[common_genes,]
      temp[temp > 0] <- 1
      mega.syn[i,] <- colSums(temp)

      temp <- input.data$non.syn[common_genes,]
      temp[temp > 0] <- 1
      mega.non.syn[i,] <- colSums(temp)
    }
  }

  mega.syn <- as.data.frame(mega.syn)
  colnames(mega.syn) <- common_sample_list
  rownames(mega.syn) <- mega.gene.names

  mega.non.syn <- as.data.frame(mega.non.syn)
  colnames(mega.non.syn) <- common_sample_list
  rownames(mega.non.syn) <- mega.gene.names
  return(list(mega.syn = mega.syn, mega.non.syn = mega.non.syn, missing.genes = missing_genes))
}

burden.test <- function(sample.list, gene.set, syn.mat, non.syn.mat, X = NULL, Y = NULL) {
  if (is.null(Y)) {
    if (is.element("Status",names(sample.list))) {
      Y = Y <- sample.list[,'Status'] - 1
    } else {
      stop("default Y variable was not found")
    }
  }
  X0 <- NULL
  if (!is.null(X)) {
    if (length(setdiff(X, names(sample.list))) > 0) {
      stop("Not all covariates are found in sample list file")
    }
    X0 <- sample.list[,X]
  }
  data.base <- cbind(Y,X0)

  p <- dim(data.base)[2] + 2
  mega.gene.names <- lapply(gene.set, function (x) x[[1]])

  p_values <- matrix(NA,length(gene.set),p)
  rownames(p_values) <- mega.gene.names
  colnames(p_values) <- c('(Intercept)', 'non.syn', 'syn', names(data.base)[-1])

  fm <- as.formula(paste("Y~", paste(colnames(p_values)[-1], collapse="+")))
  for (i in 1:length(gene.sets)) {
    if (i %% 10 == 0) {
      print(paste("Processing ", i, " out of ", length(gene.sets), " megagenes..."))
    }
    non.syn <- t(non.syn.mat[i,])
    syn <- t(syn.mat[i,])
    data.i <- cbind(data.base, non.syn , syn)
    names(data.i) <- c(names(data.base), "non.syn", "syn")
    try({
      coef <- summary(glm( fm, data = data.i, family = binomial(link = logit)))$coef
      p_values[i,rownames(coef)] <- coef[,4]
    }, silent=TRUE)
  }
  return(p_values)
}
