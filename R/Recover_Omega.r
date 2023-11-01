
###################### Step 1: Recover Omega ###################################

# refer to package lrpsadmm/R/fit_lrps.R with modification

.updatePhi.R = function(Omega, V, rho, delta){
  ## ---------------------------------------------------------------------------
  ## The name of the function: .updatePhi.R
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to update Phi in inner ADMM algorithm.
  ## ---------------------------------------------------------------------------
  Term1 = Omega - (V / rho)
  eig_Term1 = eigen(Term1, symmetric = TRUE)
  eignVal = eig_Term1$values
  eignVal[eignVal < delta] = delta
  Phi = eig_Term1$vectors %*% diag(eignVal) %*% t(eig_Term1$vectors)

  return(Phi)
}


.updateOmega.R = function(Theta, L, Phi, U, V, lambda1, mu, rho){
  ## ---------------------------------------------------------------------------
  ## The name of the function: .updateOmega.R
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to update Omega in inner ADMM algorithm.
  ## ---------------------------------------------------------------------------
  R2 = Theta - L + (U / mu)
  Term1 = R2 + V + rho * Phi
  Term2 = abs(Term1) - (lambda1 / mu)
  Term2[Term2 < 0] = 0
  Term2 = Term2 * sign(Term1)
  diag(Term2) = diag(Term1)
  Omega = (1 / (1 + rho)) * Term2

  return(Omega)
}


.updateV.R = function(Phi, Omega, V, rho){
  ## ---------------------------------------------------------------------------
  ## The name of the function: .updateV.R
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to update V in inner ADMM algorithm.
  ## ---------------------------------------------------------------------------
  V = V + rho * (Phi - Omega)

  return(V)
}

.obj_func2.R = function(Theta, L, U, Omega, lambda1, mu){
  ## ---------------------------------------------------------------------------
  ## The name of the function: .obj_func2.R
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function is corresponding to the objective function in inner ADMM
  ## algorithm.
  ## ---------------------------------------------------------------------------
  R2 = Theta - L + (U / mu)
  objval = (lambda1 / mu) * (sum(abs(Omega)) - sum(abs(diag(Omega)))) + (1 / 2) * (norm(Omega-R2, 'F'))^2

  return(objval)
}

################################################################################

.updateTheta.R = function(Sigma, Omega, L, U, mu){
  ## ---------------------------------------------------------------------------
  ## The name of the function: .updateTheta.R
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to update Theta in outer ADMM algorithm.
  ## ---------------------------------------------------------------------------
  R1 = mu * (Omega + L) - Sigma - U
  Term1 = R1 %*% R1 + 4 * mu * diag(dim(R1)[1])
  eig_Term1 = eigen(Term1, symmetric = TRUE)
  sqrt_Term1 = eig_Term1$vectors %*% diag(sqrt(eig_Term1$values)) %*% t(eig_Term1$vectors)
  Theta = (R1 + sqrt_Term1) / (2 * mu)

  return(Theta)
}


.updateOmega.admm.R = function(Theta, L, U, lambda1, mu, rho, delta, init, maxiter, abs_tol, rel_tol,
                               print_progress, print_every){
  ## ---------------------------------------------------------------------------
  ## The name of the function: .updateOmega.admm.R
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to update Omega via ADMM algorithm (inner ADMM algorithm).
  ## ---------------------------------------------------------------------------
  p = dim(Theta)[1]
  if (is.null(init)) {

    Phi = diag(p)
    Omega = diag(p)
    V = Phi * 0.0
  }else{
    parameters = init
    Phi = init$Phi
    Omega = init$Omega
    V = init$V
  }

  history = matrix(NA, nrow=0, ncol=6)
  colnames(history) = c('Iteration', 'Objval', 's_norm', 'r_norm', 'eps_pri', 'eps_dual')
  parameters = list()
  parameters$termcode = -1
  parameters$termmsg = "Maximum number of iterations reached."
  objval = .obj_func2.R(Theta, L, U, Omega, lambda1, mu)
  for (i in 1:maxiter){

    # Update Phi
    Phi = .updatePhi.R(Omega, V, rho, delta)

    # Update Omega
    Omega_old = Omega
    Omega = .updateOmega.R(Theta, L, Phi, U, V, lambda1, mu, rho)

    # Update V
    V = .updateV.R(Phi, Omega, V, rho)

    # Diagnostics
    objval_old = objval
    objval = .obj_func2.R(Theta, L, U, Omega, lambda1, mu)

    r_norm = norm(Phi - Omega, 'F')
    s_norm = norm(rho*(Omega - Omega_old), 'F')
    eps_pri = p * abs_tol + rel_tol * max(norm(Phi, 'F'), norm(Omega, 'F'))
    eps_dual = p * abs_tol + rel_tol * norm(rho * V, 'F')

    history = rbind(history, c(i, objval, s_norm, r_norm, eps_pri, eps_dual))


    if ( (r_norm < eps_pri) && (s_norm < eps_dual)){
      parameters$termcode = 0
      parameters$termmsg = 'Convergence Reached.'
      break()
    }


    if ( (print_progress) & (i > print_every)){
      if ( (i %% print_every) == 0){
        print(paste(c("Iteration:", "Obj. fun.:", "s_norm:", "r_norm:", "eps_pri:", "eps_dual:"), history[i,]))
      }
    }


  }

  parameters$Phi = Phi
  parameters$Omega = Omega
  parameters$history = as.data.frame(history)

  return(parameters)
}



.updateL.R = function(Theta, Omega, U, lambda2, mu){
  ## ---------------------------------------------------------------------------
  ## The name of the function: .updateL.R
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to update L in outer ADMM algorithm.
  ## ---------------------------------------------------------------------------
  Term1 = Theta - Omega + (U / mu)
  SVD_Term1 = svd(Term1)
  singVal = SVD_Term1$d - (lambda2 / mu)
  singVal[singVal < 0] = 0
  L = SVD_Term1$u %*% diag(singVal) %*% t(SVD_Term1$v)

  return(L)
}


.updateU.R = function(Theta, Omega, L, U, mu){
  ## ---------------------------------------------------------------------------
  ## The name of the function: .updateU.R
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to update U in outer ADMM algorithm.
  ## ---------------------------------------------------------------------------
  U = U + mu * (Theta - Omega - L)

  return(U)
}


.obj_func.R = function(Sigma, Theta, Omega, L, lambda1, lambda2){
  ## ---------------------------------------------------------------------------
  ## The name of the function: .obj_func.R
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function is corresponding to the objective function in outer ADMM
  ## algorithm.
  ## ---------------------------------------------------------------------------
  evals = eigen(Theta, symmetric = TRUE, only.values = TRUE)$values
  evals[abs(evals) < 1e-9] = 1e-8
  if (any(evals < 0)){
    return(NaN)
  }

  diag(Omega) = 0
  svd_L = svd(L)
  svd_L_values = svd_L$d

  objval = sum(diag(Sigma %*% Theta)) - sum(log(evals)) + lambda1 * sum(abs(Omega)) + lambda2 * sum(svd_L_values)

  return(objval)
}




#' Function to recover the Omega
#'
#' @usage Recover_Omega.R(Sigma, lambda1, lambda2, init, maxiter_inner, maxiter_outer, mu, rho,
#'                        delta, abs_tol, rel_tol,
#'                        print_progress_inner, print_every_inner,
#'                        print_progress_outer, print_every_outer)
#'
#' @description The function to recover the Omega via ADMM algorithm.
#'
#' @param Sigma The sample covariance matrix
#' @param lambda1 Regularization parameter corresponding to the sparsity of Omega
#' @param lambda2 Regularization parameter corresponding to the low-rank of L
#' @param init The initial value for ADMM algorithm
#' @param  maxiter_inner The maximum iteration of inner ADMM algorithm
#' @param maxiter_outer The maximum iteration of outer ADMM algorithm
#' @param mu Stepsize of the ADMM algorithm
#' @param rho Stepsize of the ADMM algorithm
#' @param delta Tolerance for updating Phi
#' @param abs_tol Absolute tolerance to stop the algorithm
#' @param rel_tol Relative tolerance to stop the algorithm
#' @param print_progress_inner  Whether the inner ADMM algorithm should report its progress
#' @param print_every_inner How often should the inner ADMM algorithm report its progress
#' @param print_progress_outer Whether the outer ADMM algorithm should report its progress
#' @param print_every_outer How often should the outer ADMM algorithm report its progress
#'
#' @return A list including all parameter matrices, the value of objective function and history
#'
#' @export
#'
Recover_Omega.R = function(Sigma, lambda1, lambda2, init, maxiter_inner, maxiter_outer, mu, rho, delta,
                           abs_tol, rel_tol,
                           print_progress_inner, print_every_inner,
                           print_progress_outer, print_every_outer){
  ## ---------------------------------------------------------------------------
  ## The name of the function: Recover_Omega.R
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to recover the Omega.
  ## ---------------------------------------------------------------------------
  p = dim(Sigma)[1]

  if (is.null(init)){
    Omega = diag(p)
    L = Omega * 0.0
    U = Omega * 0.0
  } else{
    parameters = init
    Omega = init$Omega
    L = init$L
    U = init$U
  }

  history = matrix(NA, nrow=0, ncol=6)
  colnames(history) = c('Iteration', 'Objval', 's_norm', 'r_norm', 'eps_pri', 'eps_dual')
  parameters = list()
  parameters$termcode = 1
  parameters$termmsg = "Maximum number of iterations reached."

  for (i in 1:maxiter_outer){

    # Update Theta
    Theta = .updateTheta.R(Sigma, Omega, L, U, mu)

    # Update Omega
    Omega_old = Omega
    Omega = .updateOmega.admm.R(Theta, L, U, lambda1, mu, rho, delta, init=NULL, maxiter_inner, abs_tol, rel_tol,
                                           print_progress=print_progress_inner, print_every=print_every_inner)$Omega

    # Update L
    L_old = L
    L = .updateL.R(Theta, Omega, U, lambda2, mu)

    # Update U
    U = .updateU.R(Theta, Omega, L, U, mu)

    # Diagnostics
    objval = .obj_func.R(Sigma, Theta, Omega, L, lambda1, lambda2)

    r_norm = norm(Theta - (Omega + L), 'F')
    s_norm = norm(mu * ((Omega + L) - (Omega_old + L_old)), 'F')
    eps_pri = p * abs_tol + rel_tol * max(norm(Theta, 'F'), norm(Omega + L, 'F'))
    eps_dual = p * abs_tol + rel_tol * norm(mu * U, 'F')

    history = rbind(history, c(i, objval, s_norm, r_norm, eps_pri, eps_dual))

    if ( (r_norm < eps_pri) && (s_norm < eps_dual)){
      parameters$termcode = 0
      parameters$termmsg = 'Convergence Reached.'
      break()
    }
    if ((print_progress_outer) & (i > print_every_outer)){
      if ((i %% print_every_outer) == 0){
        print(paste(c("Iteration:", "Obj. fun.:", "s_norm:", "r_norm:", "eps_pri:", "eps_dual:"), history[i,]))
      }
    }

  }

  parameters$Theta = Theta
  parameters$Omega = Omega
  parameters$L = L

  parameters$U = U

  parameters$objval = objval

  parameters$history = as.data.frame(history)

  return(parameters)
}



