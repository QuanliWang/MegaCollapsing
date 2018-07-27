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
