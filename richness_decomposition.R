library(tidyverse)
bbs <- read_csv("bbs_abundances_by_site.csv")

bbs %>%
  select(site, species, count) %>%
  spread(key = species, value = count) -> species_by_site

#Species cooccurence matrix using sorensen-based dissimilarity coefficient
prac <- matrix(rbinom(24, 1, .5), ncol = 4)
pracdf <- as.data.frame(prac)
names(pracdf) <- c('a','b','c','d')

a_matrix <- matrix(nrow = ncol(pracdf), ncol = ncol(pracdf))
bc_matrix <- matrix(nrow = ncol(pracdf), ncol = ncol(pracdf))

for (i in 1:ncol(pracdf)){
  df <- subset(pracdf, pracdf[,i] == 1)
  a = colSums(df[,i]==df)
  a_matrix[i,] <- a
  b = dim(df)[1]-a
  bc_matrix[i,] <- b
}

dissimilarity <- matrix(nrow = ncol(pracdf), ncol = ncol(pracdf))

for (i in 1:ncol(pracdf)){
  for(j in 1:ncol(pracdf)){
    b_c <- bc_matrix[i,j] + bc_matrix[j,i]
    dissimilarity[i,j] <- b_c
    dissimilarity[j,i] <- b_c
  }
}

sorenson <- dissimilarity/((2*a_matrix)+dissimilarity)

