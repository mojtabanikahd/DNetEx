library("GGMselect")
library("igraph")
library("DNetFinder")


# create two complete precision matrix
generate_reference_models <- function(number_of_nodes, number_of_samples, number_of_changes,
                                      type="ScaleFree", density_of_graph = 0.2,
                                      power = 1, mult=1)
{
  #######################################
  ########### Generate Model ############
  #######################################
  
  # zero matrix generation
  A <- matrix(rep(0,number_of_nodes*number_of_nodes), nrow=number_of_nodes)
  B <- matrix(rep(0,number_of_nodes*number_of_nodes), nrow=number_of_nodes)
  total_possible_edges <- choose(number_of_nodes, 2)
  
  # create common model
  random_weights <- rnorm(total_possible_edges)
  random_weights <- random_weights + mult*sign(random_weights)
  
  
  # define base of precision matrix
  A[lower.tri(A)] <- random_weights
  A <- t(A)
  A[lower.tri(A)] <- random_weights
  
  B <- A
  
  
  change_mask <- matrix(data = 0, nrow = number_of_nodes,
                        ncol = number_of_nodes)
  indices <- 1:choose(number_of_nodes,2)
  change_mask[lower.tri(change_mask)] <- indices
  change_mask <- t(change_mask)
  change_mask[lower.tri(change_mask)] <- indices
  mask <- sample(indices,number_of_changes)
  change_mask[change_mask %in% mask] <- -1
  change_mask[change_mask != -1] <- 0
  # create structure
  if(type == "Full") {
    A <- A*matrix(1, nrow = number_of_nodes, ncol = number_of_nodes)
    B <- B*(matrix(1, nrow = number_of_nodes, ncol = number_of_nodes)+change_mask)
  }  else if(type == "Erdos") {
    library(igraph)
    g <- erdos.renyi.game(number_of_nodes, density_of_graph, loops = F, directed = F)
    graphMatrix <- as_adjacency_matrix(g, sparse = F);
    A <- A*graphMatrix;
    B <- B*(graphMatrix+change_mask);
  }  else if(type == "ScaleFree") {
    g <- barabasi.game(number_of_nodes, power, directed = F);
    graph_matrix <- as_adjacency_matrix(g, sparse = F);
    A <- A*graph_matrix;
    B <- B*(graph_matrix+change_mask);
  } else {
    errorCondition(message = "The type value is not standard!")
  }
  
  
  # making positive definite
  minimum_eigen_value <- min(c(eigen(A)$values, eigen(B)$values))
  
  if(minimum_eigen_value < 1)
  {
    A <- A + diag(1-minimum_eigen_value, nrow = number_of_nodes,
                  ncol = number_of_nodes)
    B <- B + diag(1-minimum_eigen_value, nrow = number_of_nodes,
                  ncol = number_of_nodes)
  }
  
  
  #########################################
  ########### Generate Samples ############
  #########################################
  covariance_matrix_A <- solve(A)
  samples_A <- rmvnorm(number_of_samples, mean=rep(0,number_of_nodes), sigma=covariance_matrix_A)
  covariance_matrix_B <- solve(B)
  samples_B <- rmvnorm(number_of_samples, mean=rep(0,number_of_nodes), sigma=covariance_matrix_B)
  
  list(samples_A = samples_A, precision_matrix_A = A, covariance_matrix_A = covariance_matrix_A,
       samples_B = samples_B, precision_matrix_B = B, covariance_matrix_B = covariance_matrix_B)
}


DNetFinder_Liu2017 <- function(SA, SB, alphas, delta_star){
  est_coefGGM1=lassoGGM(SA)
  est_coefGGM2=lassoGGM(SB)
  results <- NULL
  for(alpha in alphas){
    est_DNGGM=DNetGGM(SA,SB,est_coefGGM1,est_coefGGM2,alpha)
    
    realDiffSupport <- delta_star[upper.tri(delta_star)]
    estimatedDiffSupport <- est_DNGGM[upper.tri(est_DNGGM)]

    realDiffSupport[realDiffSupport != 0] <- 1
    estimatedDiffSupport[estimatedDiffSupport != 0] <- 1

    # support evaluator
    NT <- sum(abs(realDiffSupport) == abs(estimatedDiffSupport))
    NTN <- sum(abs(realDiffSupport) == abs(estimatedDiffSupport) & realDiffSupport == 0)
    NTP <- NT - NTN
    
    NF <- sum(abs(realDiffSupport) != abs(estimatedDiffSupport))
    NFP <- sum(abs(realDiffSupport) != abs(estimatedDiffSupport) & realDiffSupport == 0)
    NFN <- NF - NFP
    
    NTPR <- NTP/(NTP+NFN)
    NFPR <- NFP/(NTN+NFP)
    
    NACC <- NT/(NT+NF)
    
    NPrecision <- NA
    if(NTP+NFP > 0)
      NPrecision <- NTP/(NTP+NFP)
    NRecall <- NTP/(NTP+NFN)
    
    FPR = NFP/(NFP+NTN)
    TPR = NTP/(NTP+NFN)
    FDR = 0
    if(NFP+NTP > 0){
      FDR = NFP / (NFP+NTP)
    }
    
    results <- rbind(results, data.frame(FPR=FPR, TPR=TPR, FDR=FDR, NPrecision=NPrecision, NRecall=NRecall, alpha=alpha,
               NTP=NTP, NTN=NTN, NFP=NFP, NFN=NFN, NTPR=NTPR, NFPR=NFPR, NACC=NACC))
  }
  results
}


evaluation <- function(estimatedDiffSupport, realDiffSupport, alpha) {
  realDiffSupport = abs(sign(realDiffSupport))
  estimatedDiffSupport = abs(sign(estimatedDiffSupport))
  
  realDiffSupport <- realDiffSupport[upper.tri(realDiffSupport)]
  estimatedDiffSupport <- estimatedDiffSupport[upper.tri(estimatedDiffSupport)]
  
  # support evaluator
  NT <- sum(abs(realDiffSupport) == abs(estimatedDiffSupport))
  NTN <- sum(abs(realDiffSupport) == abs(estimatedDiffSupport) & realDiffSupport == 0)
  NTP <- NT - NTN
  
  NF <- sum(abs(realDiffSupport) != abs(estimatedDiffSupport))
  NFP <- sum(abs(realDiffSupport) != abs(estimatedDiffSupport) & realDiffSupport == 0)
  NFN <- NF - NFP
  
  NTPR <- NTP/(NTP+NFN)
  NFPR <- NFP/(NTN+NFP)
  
  NACC <- NT/(NT+NF)
  
  NPrecision <- NA
  if(NTP+NFP > 0)
    NPrecision <- NTP/(NTP+NFP)
  NRecall <- NTP/(NTP+NFN)
  
  FPR = NFP/(NFP+NTN)
  TPR = NTP/(NTP+NFN)
  FDR = 0
  if(NFP+NTP > 0){
    FDR = NFP / (NFP+NTP)
  }
  
  results <- data.frame(FPR=FPR, TPR=TPR, FDR=FDR, NPrecision=NPrecision, NRecall=NRecall, alpha=alpha,
                        NTP=NTP, NTN=NTN, NFP=NFP, NFN=NFN, NTPR=NTPR, NFPR=NFPR, NACC=NACC)
  results
}

# DiffNetFDR_Liu2017 <- function(SA, SB, alphas, delta_star){
#   library(DiffNetFDR)
#   X <- rbind(SA, SB)
#   n_A = nrow(SA)
#   n_B = nrow(SB)
#   group <- c(rep("A", n_A), rep("B", n_B))
#   results <- NULL
#   for(alpha in alphas){                 
#         pcor.test_Liu2017 <- DiffNet.FDR(X, group, alpha, "pcor")
#         results <- rbind(results, evaluation(pcor.test_Liu2017$Diff.edge, delta_star, alpha))
#   }
#   results
# }

DiffNetFDR_Xia2015 <- function(SA, SB, alphas, delta_star){
  library(DiffNetFDR)
  X <- rbind(SA, SB)
  n_A = nrow(SA)
  n_B = nrow(SB)
  group <- c(rep("A", n_A), rep("B", n_B))
  results <- NULL
  pmat.test_Xia2015 <- DiffNet.FDR(X, group, alphas, "pmat")
  for(i in 1:length(alphas)){                  
    results <- rbind(results, evaluation(pmat.test_Xia2015[[i]]$Diff.edge, delta_star, alphas[i]))
  }
  results
}
