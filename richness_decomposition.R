library(tidyverse)
bbs <- read_csv("bbs_abundances_by_site.csv")

species_by_site <- as.data.frame(matrix(0, ncol = length(unique(bbs$species)), nrow = length(unique(bbs$site))))
names(species_by_site) <- as.integer(unique(bbs$species))

#Species by site matrix
site_ID <- unique(bbs$site)
for (i in 1:length(site_ID)) {
  site_data <- subset(bbs, site == site_ID[i])
  for (j in 1:length(site_data)){
    species <- toString(site_data$species[j])
    species_by_site[i,species] <- 1
  }
}

# bbs %>%
#   for (i in 1:3){
#     filter(bbs$site == site_ID[i]) %>%
#       for (j in 1:length(site)){
#         species <- toString(site[j])
#         species_by_site$species[i] <- 1
#       }
#   }

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

sorenson <- dissimilarity/((2*a_matrix)+dissimilarity))

