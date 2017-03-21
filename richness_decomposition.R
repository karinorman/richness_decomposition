library(tidyverse)
library(EcoSimR)

bbs <- read_csv("bbs_abundances_by_site.csv")

bbs %>%
  select(site, species, count) %>%
  spread(key = species, value = count) %>%
  select(-site) -> species_by_site
  
#replace abundance with occurence
bbs_cooc <- as_data_frame(matrix(0, nrow = nrow(species_by_site), ncol = ncol(species_by_site)))
for (i in 1:ncol(bbs_cooc)){
  locs <- which(!is.na(species_by_site[,i]))
  for (j in 1:length(locs)){
    bbs_cooc[j,i] <- 1
  }
}

#Species cooccurence matrix using sorensen-based dissimilarity coefficient
get_sorenson_matrix <- function(cooccurrence_df) {
  a_matrix <-
    matrix(nrow = ncol(cooccurrence_df),
           ncol = ncol(cooccurrence_df))
  bc_matrix <-
    matrix(nrow = ncol(cooccurrence_df),
           ncol = ncol(cooccurrence_df))
  
  for (i in 1:ncol(cooccurrence_df)) {
    df <- subset(cooccurrence_df, cooccurrence_df[, i] == 1)
    a = colSums(df[, i] == df)
    a_matrix[i, ] <- a
    b = dim(df)[1] - a
    bc_matrix[i, ] <- b
  }
  
  dissimilarity <-
    matrix(nrow = ncol(cooccurrence_df),
           ncol = ncol(cooccurrence_df))
  
  for (i in 1:ncol(cooccurrence_df)) {
    for (j in 1:ncol(cooccurrence_df)) {
      b_c <- bc_matrix[i, j] + bc_matrix[j, i]
      dissimilarity[i, j] <- b_c
      dissimilarity[j, i] <- b_c
    }
  }
  
  sorenson <- dissimilarity / ((2 * a_matrix) + dissimilarity)
  return(as_data_frame(sorenson))
}

#Test Case
prac <- matrix(rbinom(24, 1, .5), ncol = 4)
pracdf <- as.data.frame(prac)
names(pracdf) <- c('a','b','c','d')
prac_sor <- get_sorenson_matrix(pracdf)
prac_null <- cooc_null_model(t(pracdf), algo = "sim9", nReps = 1000, saveSeed = TRUE)
plot(prac_null)

#BBS Data
bbs_sor <- get_sorenson_matrix(bbs_cooc)
