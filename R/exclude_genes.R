get.pheno.geno.map <- function(files,threshold.p) {
  exclude <- NULL
  for (f in files) {
    current <- fread(f, header=TRUE, stringsAsFactors = FALSE)
    exclude <- rbind(exclude,current)
  }
  exclude <- exclude %>% filter(pValue <= threshold.p)
  exclude <- exclude %>% dplyr::select(c("phenotype","genotype"))
  exclude <- exclude %>% group_by(phenotype) %>% mutate(genolist = paste(genotype, collapse = ",")) %>%
    dplyr::select(phenotype,genolist) %>%unique()
  return(exclude)
}
