library("GGMselect")
library("igraph")
library("DNetFinder")

DNetFinder_Liu2017_method <- function(SA, SB, alpha){
  library(DNetFinder)
  est_coefGGM1=lassoGGM(SA)
  est_coefGGM2=lassoGGM(SB)
  est_DNGGM=DNetGGM(SA,SB,est_coefGGM1,est_coefGGM2,alpha)
  est_DNGGM
}

DiffNetFDR_Liu2017_method <- function(SA, SB, alpha){
  library(DiffNetFDR)
  X <- rbind(SA, SB)
  n_A = nrow(SA)
  n_B = nrow(SB)
  group <- c(rep("A", n_A), rep("B", n_B))
  res  <- tryCatch(
    {
      pcor.test_Liu2017 <- DiffNet.FDR(X, group, alpha, "pcor")
      pcor.test_Liu2017[[1]]
    },
    error = function(e) {
      print("No differential edges identified.")
      print(substr(conditionMessage(e), 1, 100))
      matrix(0, ncol(SA), ncol(SA))
    }
  )
  res
}

DiffNetFDR_Xia2015_method <- function(SA, SB, alpha){
  library(DiffNetFDR)
  X <- rbind(SA, SB)
  n_A = nrow(SA)
  n_B = nrow(SB)
  group <- c(rep("A", n_A), rep("B", n_B))
  res  <- tryCatch(
    {
      pmat.test_Xia2015 <- DiffNet.FDR(X, group, alpha, "pmat")
      pmat.test_Xia2015[[1]]
    },
    error = function(e) {
      print("No differential edges identified.")
      print(substr(conditionMessage(e), 1, 100))
      matrix(0, ncol(SA), ncol(SA))
    }
  )
  res
}