exclude.genes <- function(input.data, exclude.file) {
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

as.mega.matrix<- function(gene.sets, input.data) {
  mega.gene.names<- sapply(gene.sets, function(x) x[1])
  common_sample_list <- rownames(input.data$sample.list)
  has.syn <- "syn" %in% names(input.data)
  if (has.syn) {
    mega.syn <- matrix(NA,nrow = length(gene.sets), ncol = length(common_sample_list))
  }
  mega.non.syn <- matrix(NA,nrow = length(gene.sets), ncol = length(common_sample_list))
  missing_genes = list()
  for (i in 1:length(gene.sets)) {
    if (i %% 10 == 0) {
      print(paste("Processing ", i, " out of ", length(gene.sets), " megagenes..."))
    }
    genes <- gene.sets[[i]][3:length(gene.sets[[i]])]
    common_genes <- intersect(genes,rownames(input.data$non.syn))
    missing_genes[[mega.gene.names[i]]] <- setdiff(genes,rownames(input.data$non.syn))
    if (length(common_genes)>0) {
      if (has.syn) {
        temp <- input.data$syn[common_genes,]
        temp[temp > 0] <- 1
        mega.syn[i,] <- colSums(temp)
      }
      temp <- input.data$non.syn[common_genes,]
      temp[temp > 0] <- 1
      mega.non.syn[i,] <- colSums(temp)
    }
  }
  if (has.syn) {
    mega.syn <- as.data.frame(mega.syn)
    colnames(mega.syn) <- common_sample_list
    rownames(mega.syn) <- mega.gene.names
  }
  mega.non.syn <- as.data.frame(mega.non.syn)
  colnames(mega.non.syn) <- common_sample_list
  rownames(mega.non.syn) <- mega.gene.names

  adjusted.gene.set.size <- unlist(lapply(gene.sets,length))  - 2 - unlist(lapply(missing_genes,length))

  if (has.syn) {
    return(list(mega.syn = mega.syn, mega.non.syn = mega.non.syn, missing.genes = missing_genes,
                adjusted.gene.set.size = adjusted.gene.set.size))
  } else {
    return(list(mega.non.syn = mega.non.syn, missing.genes = missing_genes,
                adjusted.gene.set.size = adjusted.gene.set.size))
  }
}

burden.test <- function(sample.list, gene.sets, syn.mat, non.syn.mat, X = NULL, Y = NULL) {
  if (is.null(Y)) {
    if (is.element("Status",names(sample.list))) {
      Y <- sample.list[,'Status'] - 1
    } else {
      stop("default Y variable was not found")
    }
  }
  median.syn.case <- median(rowSums(syn.mat[,Y == 1]))
  median.syn.control <- median(rowSums(syn.mat[,Y == 0]))
  median.non.syn.case <- median(rowSums(non.syn.mat[,Y == 1]))
  median.non.syn.control <- median(rowSums(non.syn.mat[,Y == 0]))


  X0 <- NULL
  if (!is.null(X)) {
    if (length(setdiff(X, names(sample.list))) > 0) {
      stop("Not all covariates are found in sample list file")
    }
    X0 <- sample.list[,X]
  }
  data.base <- cbind(Y,X0)
  if (length(X) == 1) {
    data.base <- as.data.frame(data.base)
    names(data.base) <- c("Y",X)
  }

  p <- dim(data.base)[2] + 2
  mega.gene.names <- lapply(gene.sets, function (x) x[[1]])

  p_values <- matrix(NA,length(gene.sets),p)
  rownames(p_values) <- mega.gene.names
  colnames(p_values) <- c('(Intercept)', 'non.syn', 'syn', names(data.base)[-1])

  betas <- matrix(NA,length(gene.sets),p)
  rownames(betas) <- mega.gene.names
  colnames(betas) <- c('(Intercept)', 'non.syn', 'syn', names(data.base)[-1])

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
      betas[i,rownames(coef)] <- coef[,1]
    }, silent=TRUE)
  }
  return(list(p_values = p_values, betas = betas, median.syn.case = median.syn.case,
              median.syn.control = median.syn.control, median.non.syn.case = median.non.syn.case,
              median.non.syn.control = median.non.syn.control))
}
