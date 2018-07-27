get.pvalues_logit <-function(Y, V1, gender, sV1, n.permutations = 1000) {
  #Number of cases and controls
  n.samples <- length(Y)
  n.cases <- sum(Y==1)
  n.controls <- n.samples - n.cases

  Y.perm.template = numeric(n.samples)

  #permutation, save all p-values just in case median will be needed later on
  p.observed <- 1.0 #assume not significant
  try({
    p.observed <- summary(glm( Y ~ V1 + gender + sV1, family = binomial(link = logit)))$coef["V1",4]
  }, silent=TRUE)

  P.Values <- numeric(n.permutations) + 1
  for (i in 1: n.permutations) {
    Y.perm <- Y.perm.template
    Y.perm[sample.int(n.samples, size = n.cases, replace = FALSE)] <-1
    try({
      P.Values[i] <- summary(glm( Y.perm ~ V1 + gender + sV1, family = binomial(link = logit)))$coef["V1",4]
    }, silent=TRUE)
  }
  P.perm <- length(which(P.Values<p.observed)) / n.permutations

  out <- list()
  out$perm <- P.perm
  out$observed <- p.observed
  out
}

get.pvalues_logit_null <-function(Y, gender,non_syn,syn, n.permutations = 100) {

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
  out <- list()
  out$perm <- P.perm
  out$observed <- p.observed
  out$matrix = permed_p_values
  out$perm.labels <- KK
  out
}
