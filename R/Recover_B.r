#' Function to recover B
#'
#' @usage Recover_B.R(X, LayerOrder, kappa_n, nu_n)
#'
#' @description The function to recover/estimate B.
#'
#' @param X The data matrix, n*p matrix
#' @param  LayerOrder The causal ordering among chain components
#' @param kappa_n A hard threshold for truncated singular value decomposition(SVD)
#' @param  nu_n An element-wise hard threshold for estimated B
#'
#' @return The estimated B, p*p matrix
#'
#' @export
#'
#' @import MASS


################################################################################
############################ Recover B #########################################
################################################################################

Recover_B.R = function(X, LayerOrder, kappa_n, nu_n){
  ## ---------------------------------------------------------------------------
  ## The name of the function: Recover_B.R
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to recover/estimate B.
  ## ---------------------------------------------------------------------------
  p = ncol(X)
  numLayer = length(LayerOrder)

  # regression
  B_reg = matrix(0,p,p)
  Regressor = NULL

  for (i in 2:numLayer){
    Response = LayerOrder[[i]]
    Regressor = union(Regressor,LayerOrder[[i-1]])
    B_reg[Response,Regressor] = MASS::ginv(t(X[,Regressor])%*%X[,Regressor])%*%t(X[,Regressor])%*%X[,Response]
  }
  zero_state = which(B_reg==0)

  # sparse SVD
  svd_B_reg = svd(B_reg)
  d_B_reg = svd_B_reg$d
  d_B_reg[d_B_reg < kappa_n] = 0

  B_ssvd = svd_B_reg$u%*%diag(d_B_reg)%*%t(svd_B_reg$v)

  # element-wise hard thresholding
  B = B_ssvd
  B[zero_state] = 0
  B[abs(B) < nu_n] = 0


  return(B)
}
