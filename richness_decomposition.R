library(tidyverse)
bbs <- read_csv("bbs_abundances_by_site.csv")

species_by_site <- as.data.frame(matrix(0, ncol = length(unique(bbs$species)), nrow = length(unique(bbs$site))))
names(species_by_site) <- as.integer(unique(bbs$species))

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

