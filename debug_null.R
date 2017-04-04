############# Examples from the documentation of cooc_null_model ####################################

library(EcoSimR)
library(tidyverse)

## Example is not identical if we do not save seed
finchMod <- cooc_null_model(dataWiFinches, algo="sim9",nReps=10000,burn_in = 500)
finch2 <- cooc_null_model(dataWiFinches, algo="sim9",nReps=10000,burn_in = 500)

identical(finchMod, finch2)




## Example that is repeatable with a saved seed
finchMod <- cooc_null_model(dataWiFinches, algo="sim1",saveSeed = TRUE)
a <- mean(finchMod$Sim)
## Run the model with the seed saved

finchMod <- cooc_null_model(dataWiFinches, algo="sim1",saveSeed=T)
## Check model output
b <- mean(finchMod$Sim)

## So much for the documentation, these are still not identical.
identical(a, b)

## This doesn't even run, just throws an error
## reproduce_model(finchMod$Sim)


## Not even sure why this is included, but it is not identical
finchMod <- cooc_null_model(dataWiFinches, algo="sim1")
mean(finchMod$Sim)
## reproduce_model(finchMod$Sim)



############################ Example from Kari's code ###########################

prac <- matrix(rbinom(24, 1, .5), ncol = 4)
pracdf <- as.data.frame(prac)
names(pracdf) <- c('a','b','c','d')

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

prac_sor <- get_sorenson_matrix(pracdf)

df <- as.data.frame(prac_sor)


## These are identical:
n1 <- cooc_null_model(df, algo = "sim9", nReps = 1000, saveSeed = FALSE)$Randomized.Data
n2 <- cooc_null_model(df, algo = "sim9", nReps = 1000, saveSeed = FALSE)$Randomized.Data
identical(n1, n2)

## All identical except time stamp:
n1 <- cooc_null_model(df, algo = "sim9", nReps = 1000, saveSeed = FALSE)
n2 <- cooc_null_model(df, algo = "sim9", nReps = 1000, saveSeed = FALSE)
map2(n1, n2, identical)

## Randomized.Data is not identical with sim1:
n1 <- cooc_null_model(df, algo = "sim1", nReps = 1000, saveSeed = FALSE)$Randomized.Data
n2 <- cooc_null_model(df, algo = "sim1", nReps = 1000, saveSeed = FALSE)$Randomized.Data
identical(n1, n2)


## Okay, so things are looking pretty hokey at this point.... Examining cooc_null_model shows completely different dispatch methods for sim9 vs the rest:
cooc_null_model

## Calling sim9 routine directly, we find that these are still identical
n1 <- sim9(df, algo = "sim9", metric = "c_score")$Randomized.Data
n2 <- sim9(df, algo = "sim9", metric = "c_score")$Randomized.Data
identical(n1,n2)

## So time to dig into code for sim9
sim9

## Looks like it calls sim9_single, whatever that is.  Stipping out that part of the code, we
## see sim9_single is also giving identical values on each call:

df <- speciesData
metricF <- get("c_score")
Obs <- metricF(as.matrix(speciesData))
msim <- speciesData[rowSums(speciesData) > 0, ]
n1 <- sim9_single(msim)
n2 <- sim9_single(msim)

identical(n1, n2)

ex1 <- matrix(rbinom(100, 1, 0.5), nrow = 10)
## This is not expected, or at least it doesn't occur with the default data of the function:
identical(sim9_single(ex1), sim9_single(ex1))

## So, what's special about df?  Perhaps something in the conversions of speciesData is causing this...



