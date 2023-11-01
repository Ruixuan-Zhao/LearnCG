#' Evaluation function for the estimated chain graph
#'
#' @usage evaluation.metric.final(est_Omega, est_B, true_Omega, true_B)
#'
#' @description The function to evaluate the numerical performances of methods
#' in terms of the estimation accuracy of undirected edges, directed edges and
#' the overall chain graph.
#'
#' @param est_Omega The estimated Omega
#' @param  est_B The estimated B
#' @param true_Omega The true Omega
#' @param  true_B The true B
#'
#' @return A list of evaluation results: Recall, Precision, MCC and SHD
#'
#' @export
#'
#' @import EvaluationMeasures
#'
#' @examples
#' \donttest{
#' evaluation.metric.final(est_Omega=matrix(rbinom(100,1,0.5),10,10),
#'          est_B=diag(10), true_Omega=matrix(rbinom(100,1,0.5),10,10),
#'          true_B=diag(10))
#' }
#'






################################################################################
########################## Evaluation function #################################
################################################################################

evaluation.metric.final = function(est_Omega, est_B, true_Omega, true_B){
  ## ---------------------------------------------------------------------------
  ## The name of the function: evaluation.metric.final
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to evaluate the numerical performances of methods.
  ## ---------------------------------------------------------------------------
  ## Input:
  ## @ est_Omega: the estimated Omega
  ## @ est_B: the estimated B
  ## @ true_Omega: the true Omega
  ## @ true_B: the true B
  ## ---------------------------------------------------------------------------
  ## Output:
  ## @ Eval_res_both: the list of evaluation results: Recall, Precision, MCC and SHD
  ## ---------------------------------------------------------------------------

  est_Omega_01 = est_Omega
  est_Omega_01[which(est_Omega_01 != 0)] = 1

  est_B_01 = est_B
  est_B_01[which(est_B_01 != 0)] = 1

  true_Omega_01 = true_Omega
  true_Omega_01[which(true_Omega_01 != 0)] = 1

  true_B_01 = true_B
  true_B_01[which(true_B_01 != 0)] = 1

  est_both = est_Omega + est_B
  est_both_01 = est_both
  est_both_01[which(est_both_01 != 0)] = 1

  true_both = true_Omega + true_B
  true_both_01 = true_both
  true_both_01[which(true_both_01 != 0)] = 1


  ### evaluate Omega (undirected edges)
  true_vec_UE = as.vector(true_Omega_01[which(lower.tri(true_Omega_01))])
  est_vec_UE = as.vector(est_Omega_01[which(lower.tri(est_Omega_01))])

  Recall_UE = EvaluationMeasures::EvaluationMeasures.Recall(true_vec_UE,est_vec_UE)
  Precision_UE = EvaluationMeasures::EvaluationMeasures.Precision(true_vec_UE,est_vec_UE)
  if (is.nan(Precision_UE)){Precision_UE=0}

  MCC_UE = EvaluationMeasures::EvaluationMeasures.MCC(true_vec_UE,est_vec_UE)
  if (is.nan(MCC_UE)){MCC_UE=0}

  ### evaluate B (directed edges)
  Recall_B = EvaluationMeasures::EvaluationMeasures.Recall(as.vector(true_B_01),as.vector(est_B_01))
  Precision_B = EvaluationMeasures::EvaluationMeasures.Precision(as.vector(true_B_01),as.vector(est_B_01))
  if (is.nan(Precision_B)){Precision_B=0}

  MCC_B = EvaluationMeasures::EvaluationMeasures.MCC(as.vector(true_B_01),as.vector(est_B_01)) #MCC
  if (is.nan(MCC_B)){MCC_B=0}

  # structural hamming distance
  ##############################################################################
  ###### reference: comp.cgs in https://github.com/majavid/AMPCGs2019 ##########
  ##############################################################################

  ## computing structural Hamming distance
  #Structural Hamming distance is defined as the total number of operations needed to convert one
  # graph to the other. Each operation must be one of the following: (1) add or delete an
  # edge, or (2) add, remove or reverse an orientation of an edge.
  p=nrow(est_both_01)
  shd = 0
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      if((true_both_01[i,j]!=0 || true_both_01[j,i]!=0)&&(est_both_01[i,j]==0)&&(est_both_01[j,i]==0)){
        ###missing edge in the estimated CG: add appropriate edge
        shd = shd + 1
      }
      if((est_both_01[i,j]!=0 || est_both_01[j,i]!=0)&&(true_both_01[i,j]==0)&&(true_both_01[j,i]==0)){
        ###extra edge in the estimated CG: remove it
        shd = shd + 1
      }
      if((true_both_01[i,j]+true_both_01[j,i]==1)&&(est_both_01[i,j]+est_both_01[j,i]==2)){
        ###there is a directed edge in the true CG, but the corresponding edge in estimated CG is undirected: add proper orientation
        shd = shd + 1
      }
      if((true_both_01[i,j]+true_both_01[j,i]==1)&&(est_both_01[i,j]+est_both_01[j,i]==1)&&(true_both_01[i,j]!=est_both_01[i,j])){
        ###-->/<-- in the true CG, but <--/--> in the estimated CG, respectively: reverse the orientation
        shd = shd + 1
      }
      if((true_both_01[i,j]+true_both_01[j,i]==2) && (est_both_01[i,j]+est_both_01[j,i]==1)){
        ###there is an undirected edge in the true CG, but the corresponding edge in estimated CG is directed -->/<--: remove the orientation
        shd = shd + 1
      }
    }
  }

  Eval_res_both = c(Recall_UE, Precision_UE, MCC_UE, Recall_B, Precision_B, MCC_B, shd)
  names(Eval_res_both) = c("Recall_UE", "Precision_UE", "MCC_UE", "Recall_DE", "Precision_DE",  "MCC_DE", "SHD")

  return(Eval_res_both)

}






