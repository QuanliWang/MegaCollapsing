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

as.mega.matrix<- function(gene.sets, input.data) {
  sample.list <- colnames(input.data$non.syn)
  mega.gene.names <- rownames(gene.sets$mat)

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
      input.data$syn[input.data$syn > 0] <- 1
  }
  input.data$non.syn <- input.data$non.syn[gene.index,]
  input.data$non.syn[input.data$non.syn > 0] <- 1

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

burden.test <- function(sample.list, gene.sets, syn.mat, non.syn.mat, X = NULL, Y = NULL, ncores = 1) {

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
  mega.gene.names <- rownames(gene.sets$mat)

  p_values <- matrix(NA,length(mega.gene.names),p)
  rownames(p_values) <- mega.gene.names
  colnames(p_values) <- c('(Intercept)', 'non.syn', 'syn', names(data.base)[-1])

  betas <- matrix(NA,length(mega.gene.names),p)
  rownames(betas) <- mega.gene.names
  colnames(betas) <- c('(Intercept)', 'non.syn', 'syn', names(data.base)[-1])
  fm <- as.formula(paste("Y~", paste(colnames(p_values)[-1], collapse="+")))

  if (ncores >1) {
      # Set up parallel backend
      cl.create <- FALSE
      if (foreach::getDoParWorkers() > 1) {
        if (foreach::getDoParWorkers() != ncores) {
          message(paste("Parallel backend already registered and is inconsistent with ncores."))
        }
        ncores <- getDoParWorkers()
      } else {
          cl.create <- TRUE
          cl <- parallel::makeCluster(ncores, outfile = "")
          doParallel::registerDoParallel(cl)
          on.exit({
            parallel::stopCluster(cl)
            foreach::registerDoSEQ()
          })
      }

      itercount <- iterators::icount(length(mega.gene.names))
      out <- suppressWarnings(foreach::foreach(i = itercount, .packages = "MegaCollapsing", .errorhandling="pass") %dopar% {
                if (i %% 10 == 0) {
                      print(paste("Processing ", i, " out of ", length(mega.gene.names), " megagenes..."))
                }
                non.syn <- non.syn.mat[i,]
                syn <- syn.mat[i,]
                data.i <- cbind(data.base, non.syn , syn)
                names(data.i) <- c(names(data.base), "non.syn", "syn")
                try({
                    coef <- summary(glm( fm, data = data.i, family = binomial(link = logit)))$coef
                    p_value <- coef[,4]
                    beta <- coef[,1]
                   list(p_value, beta, i)
                }, silent=TRUE)
              })
      indexes <- sapply(out,function(x) x[[3]])
      for (i in 1:length(indexes)) {
        p_values[indexes[i],names(out[[i]][[2]])] <- out[[i]][[1]]
        betas[indexes[i],names(out[[i]][[2]])] <- out[[i]][[2]]
      }
  } else {
      for (i in 1:length(mega.gene.names)) {
        if (i %% 10 == 0) {
          print(paste("Processing ", i, " out of ", length(mega.gene.names), " megagenes..."))
        }
        non.syn <- non.syn.mat[i,]
        syn <- syn.mat[i,]
        data.i <- cbind(data.base, non.syn , syn)
        names(data.i) <- c(names(data.base), "non.syn", "syn")
        try({
          coef <- summary(glm( fm, data = data.i, family = binomial(link = logit)))$coef
          p_values[i,rownames(coef)] <- coef[,4]
          betas[i,rownames(coef)] <- coef[,1]
        }, silent=TRUE)
      }
  }
  return(list(p_values = p_values, betas = betas, median.syn.case = median.syn.case,
              median.syn.control = median.syn.control, median.non.syn.case = median.non.syn.case,
              median.non.syn.control = median.non.syn.control))
}
