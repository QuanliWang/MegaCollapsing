
burden.test.with.permutation <-function(sample.list, gene.set, syn.mat, non.syn.mat, X = NULL, Y = NULL, n.permutations = 100) {
  # get observed p-values
  print("computing observed p-values...")
  p.observed <- burden.test(sample.list, gene.set, syn.mat, non.syn.mat, X, Y)[,"non.syn"]
  p.observed[is.na(p.observed)] <- 1

  if (is.null(Y)) {
    if (is.element("Status",names(sample.list))) {
      Yv <- sample.list[,'Status'] - 1
    } else {
      stop("default Y variable was not found")
    }
  }

  #Number of cases and controls
  n.samples <- length(Yv)
  n.cases <- sum(Yv==1)
  n.controls <- n.samples - n.cases
  Y.perm.template = numeric(n.samples)

  #permutation, save all p-values just in case median will be needed later on
  permed_p_values <- matrix(1,length(gene.set),n.permutations)
  KK <- matrix(0,n.permutations,n.cases)
  sample.list.perm <- sample.list
  for (i in 1: n.permutations) {
    Y.perm <- Y.perm.template
    K <- sample.int(n.samples, size = n.cases, replace = FALSE)
    Y.perm[K] <- 1
    sample.list.perm[,'Status'] <- Y.perm + 1
    print(paste("Running permutaiton ", i))
    ps <- burden.test(sample.list.perm, gene.set, syn.mat, non.syn.mat, X, Y)[,"non.syn"]
    ps[is.na(ps)] <- 1
    permed_p_values[,i] <- ps
    KK[i,] <- K
  }
  P.perm  <- apply(permed_p_values, 1, median)
  rownames(permed_p_values) <- rownames(non.syn.mat)

  P.freq <- P.perm
  for (i in 1: length(P.perm)) {
    P.freq[i] <- sum(permed_p_values[i,] < p.observed[i]) / n.permutations
  }

  out <- list()
  out$perm <- P.perm
  out$observed <- p.observed
  out$matrix = permed_p_values
  out$perm.labels <- KK

  p_summary <- matrix(NA,length(P.perm),3)
  p_summary[,1] <- P.perm
  p_summary[,2] <- p.observed
  p_summary[,3] <- P.freq
  rownames(p_summary) <-  lapply(gene.set, function (x) x[[1]])
  colnames(p_summary) <- c('perm','observed','freq')
  out$p_summary <- p_summary
  out
}

.burden.test.with.permutation.old <-function(Y, gender,non_syn,syn, n.permutations = 100) {
  # get observed p-values
  print("computing observed p-values...")
  n <- dim(non_syn)[1]
  p.observed <- numeric(n) + 1
  for (i in 1:n) {
    if (i %% 100) {}
    print(paste(i, "out of ", n))
    V1 = non_syn[i,]
    sV1 = syn[i,]
    try({
      p.observed[i] <- summary(glm( Y ~ V1 + gender + sV1, family = binomial(link = logit)))$coef["V1",4]
    }, silent=TRUE)
  }

  #Number of cases and controls
  n.samples <- length(Y)
  n.cases <- sum(Y==1)
  n.controls <- n.samples - n.cases
  Y.perm.template = numeric(n.samples)

  #permutation, save all p-values just in case median will be needed later on

  permed_p_values <- matrix(1,dim(non_syn)[1],n.permutations)
  KK <- matrix(0,n.permutations,n.cases)
  for (i in 1: n.permutations) {
    Y.perm <- Y.perm.template
    K <- sample.int(n.samples, size = n.cases, replace = FALSE)
    Y.perm[K] <-1
    for (j in 1:dim(non_syn)[1]) {
      print(c(i,j))
      try({
        V1 = non_syn[j,]
        sV1 = syn[j,]
        permed_p_values[j,i] <- summary(glm(Y.perm ~ V1 + gender + sV1, family = binomial(link = logit)))$coef["V1",4]
      }, silent=TRUE)
    }
    KK[i,] <- K
  }
  P.perm  <- apply(permed_p_values, 1, median)
  P.freq <- P.perm
  for (i in 1: length(P.perm)) {
    P.freq[i] <- sum(permed_p_values[i,] < p.observed[i]) / n.permutations
  }
  out <- list()
  out$perm <- P.perm
  out$observed <- p.observed
  out$matrix = permed_p_values
  out$perm.labels <- KK

  p_summary <- matrix(NA,length(P.perm),3)
  p_summary[,1] <- P.perm
  p_summary[,2] <- p.observed
  p_summary[,3] <- P.freq
  rownames(p_summary) <-  rownames(non_syn)
  colnames(p_summary) <- c('perm','observed', 'freq')
  out$p_summary <- p_summary
  out
}
