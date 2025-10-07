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