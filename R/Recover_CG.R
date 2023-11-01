#' Function to reconstruct the chain graph
#'
#' @usage Recover.CG(X, lambda1_n, lambda2_n, kappa_n, nu_n)
#'
#' @description The function to reconstruct the chain graph / estimate (Omega, B).
#'
#' @param X The data matrix, n*p matrix
#' @param  lambda1_n Regularization parameter corresponding to the sparsity of Omega
#' @param lambda2_n Regularization parameter corresponding to the low-rank of L
#' @param kappa_n A hard threshold for truncated singular value decomposition(SVD)
#' @param  nu_n An element-wise hard threshold for estimated B
#'
#' @return The estimated Omega and B
#'
#' @export
#'
#' @importFrom stats cov
#'
#' @examples
#'gen_data = simul.CG.Ex1(p.node = 50, n.sample = 500, p.un = 0.03,
#'                        p.hub = 0.2, p.hub.dir = 0.8, CG.seed = 1)
#'X = gen_data$X
#'est_res = Recover.CG(X = X, lambda1_n = 0.4*500^(-3/8),
#'                     lambda2_n = 0.8*500^(-3/8), kappa_n = 5.5*500^(-3/8),
#'                     nu_n = 3.8*500^(-1/4))
#'est_Omega = est_res$est_Omega
#'est_B = est_res$est_B
#'

Recover.CG = function(X, lambda1_n, lambda2_n, kappa_n, nu_n){
  ## ---------------------------------------------------------------------------
  ## The name of the function: Recover.CG
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to reconstruct the chain graph / estimate (Omega, B).
  ## ---------------------------------------------------------------------------
  ## Input:
  ## @ X: sample matrix
  ## @ lambda1_n, lambda2_n: two regularization parameters in the optimization of (6)
  ## @ kappa_n: parameter for truncated SVD
  ## @ nu_n: element-wise hard thresholding for estimated B
  ## ---------------------------------------------------------------------------
  ## Output:
  ## @ est_Omega: estimated precision matrix of noise corresponding to the undirected edges
  ## @ est_B: estimated coefficient matrix corresponding to the directed edges
  ## ---------------------------------------------------------------------------

  Sigma = cov(X)
  res_Omega = Recover_Omega.R(Sigma = Sigma, lambda1 = lambda1_n, lambda2 = lambda2_n, init = NULL, maxiter_inner = 10000, maxiter_outer = 10000, mu = 1, rho = 0.5, delta = 1e-04,
                              abs_tol = 5e-04, rel_tol = 5e-02,
                              print_progress_inner = FALSE, print_every_inner = 1,
                              print_progress_outer = FALSE , print_every_outer = 50)

  est_Omega = res_Omega$Omega

  est_LayerOrder = Recover_order.R(X = X,Omega = est_Omega)$layerOrder

  est_B = Recover_B.R(X = X, LayerOrder = est_LayerOrder, kappa_n = kappa_n, nu_n = nu_n)

  return(list(est_Omega = est_Omega, est_B = est_B))
}
