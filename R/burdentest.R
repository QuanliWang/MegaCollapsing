do.fastLR <- function(data.i, Y.index = 1) {
    YY <- data.i[,Y.index]
    XX <- as.matrix(data.i[,-Y.index])
    index <- which(is.na(YY))
    if (length(index)>0) {
      YY <- YY[-index]
      XX <- XX[-index,]
    }
    XX <- cbind(1,XX)
    res2 <- fastLR(XX, YY)
    v <- res2$fitted.values * (1 - res2$fitted.values)
    b <- diag(solve(sweep(t(XX), 2, v, "*") %*% XX))
    se_beta <- sqrt(b)
    z <- res2$coefficients / se_beta
    p <- pnorm(abs(z), lower.tail = FALSE) * 2
    return (list(p = p, beta = res2$coefficients))
}

burden.test <- function(sample.list, syn.mat, non.syn.mat, X = NULL, Y = NULL, ncores = 1, excluded_samples = NULL, do.fast = FALSE) {

  if (is.null(Y)) {
    if (is.element("Status",names(sample.list))) {
      Y <- sample.list[,'Status'] - 1
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
  if (length(X) == 1) {
    data.base <- as.data.frame(data.base)
    names(data.base) <- c("Y",X)
  }

  p <- dim(data.base)[2] + 2
  mega.gene.names <- rownames(syn.mat)

  p_values <- matrix(NA,length(mega.gene.names),p)
  rownames(p_values) <- mega.gene.names
  colnames(p_values) <- c('(Intercept)', names(data.base)[-1], 'non.syn', 'syn')

  betas <- matrix(NA,length(mega.gene.names),p)
  rownames(betas) <- mega.gene.names
  colnames(betas) <- c('(Intercept)', names(data.base)[-1], 'non.syn', 'syn')
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
      cl <- parallel::makeCluster(ncores, outfile = "", setup_timeout = 300)
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
        if (!is.null(excluded_samples)) {
          common.samples <- intersect(excluded_samples,rownames(data.i))
          if (length(common.samples) > 0) {
            common.samples.index <- match(common.samples,rownames(data.i))
            data.i <- data.i[-common.samples.index,]
          }
        }
        if (do.fast) {
            r <- do.fastLR(data.i)
            p_value <- r$p
            beta <- r$beta
        } else {
            coef <- summary(glm( fm, data = data.i, family = binomial(link = logit)))$coef
            p_value <- coef[,4]
            beta <- coef[,1]
        }
        list(p_value, beta, i)
      }, silent=TRUE)
    })
    indexes <- sapply(out,function(x) x[[3]])
    for (i in 1:length(indexes)) {
      p_values[indexes[i],] <- out[[i]][[1]]
      betas[indexes[i],] <- out[[i]][[2]]
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
        if (!is.null(excluded_samples)) {
            common.samples <- intersect(excluded_samples,rownames(data.i))
            if (length(common.samples) > 0) {
                common.samples.index <- match(common.samples,rownames(data.i))
                data.i <- data.i[-common.samples.index,]
            }
        }
        if (do.fast) {
            r <- do.fastLR(data.i)
            p_values[i,] <- r$p
            betas[i,] <- r$beta
        } else {
            coef <- summary(glm( fm, data = data.i, family = binomial(link = logit)))$coef
            p_values[i,] <- coef[,4]
            betas[i,] <- coef[,1]
        }

      }, silent=TRUE)
    }
  }
  return(list(p_values = p_values, betas = betas, mega.gene.names = mega.gene.names))
}


