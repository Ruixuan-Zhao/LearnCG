#' Function to recover the causal ordering of chain components
#'
#' @usage Recover_order.R(X,Omega)
#'
#' @description The function to recover the causal ordering of chain components.
#'
#' @param X The data matrix, n*p matrix
#' @param  Omega The precision matrix of noise Omega, p*p matrix
#'
#' @return The estimated casual ordering of chain components
#'
#' @export
#'
#' @import igraph
#'
#' @importFrom stats var cov




################################################################################
##################### Recover the topological order of layers ##################
################################################################################

Recover_order.R = function(X,Omega){
  ## ---------------------------------------------------------------------------
  ## The name of the function: Recover_order.R
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to recover the causal ordering of chain components.
  ## ---------------------------------------------------------------------------

  # identify layers (chain component)
  Omega_backpack = Omega
  diag(Omega_backpack) = 0
  undirectG = igraph::graph_from_adjacency_matrix(Omega_backpack,"undirected", weighted = TRUE)
  layersINF = igraph::components(undirectG)
  layerID = layersINF$membership
  layerNUM = layersINF$no

  D1 = rep(NA,layerNUM)
  for (j in 1:layerNUM){
    nodesInLayer = which(layerID==j)
    if (length(nodesInLayer)==1){
      D1[j] = var(X[,nodesInLayer]) - solve(Omega)[nodesInLayer, nodesInLayer]
    }else{
      D1[j] = max(diag(cov(X[,nodesInLayer]) - solve(Omega)[nodesInLayer, nodesInLayer]))
    }

  }
  LayerHOrder = sort(D1,index.return = T)$ix[1]
  LayerToOrder = setdiff(c(1:layerNUM), LayerHOrder)

  layerOrder = list(which(layerID == LayerHOrder))

  index = 1

  while(length(LayerToOrder) > 1){
    D2 = rep(NA, length(LayerToOrder))
    for (j in LayerToOrder){
      nodesInLayer = which(layerID==j)
      if (length(nodesInLayer)==1){
        if (length(LayerHOrder)==1){
          D2[j] = var(X[,nodesInLayer]) - cov(X[,nodesInLayer], X[,LayerHOrder])%*%solve(var(X[,LayerHOrder]))%*%cov(X[,LayerHOrder], X[,nodesInLayer]) - solve(Omega)[nodesInLayer, nodesInLayer]
        }else{
          D2[j] = var(X[,nodesInLayer]) - cov(X[,nodesInLayer], X[,LayerHOrder])%*%solve(cov(X[,LayerHOrder]))%*%cov(X[,LayerHOrder], X[,nodesInLayer]) - solve(Omega)[nodesInLayer, nodesInLayer]
        }

      }else{
        if (length(LayerHOrder)==1){
          D2[j] = max(diag(var(X[,nodesInLayer]) - cov(X[,nodesInLayer], X[,LayerHOrder])%*%solve(var(X[,LayerHOrder]))%*%cov(X[,LayerHOrder], X[,nodesInLayer]) - solve(Omega)[nodesInLayer, nodesInLayer]))
        }else{
          D2[j] = max(diag(var(X[,nodesInLayer]) - cov(X[,nodesInLayer], X[,LayerHOrder])%*%solve(cov(X[,LayerHOrder]))%*%cov(X[,LayerHOrder], X[,nodesInLayer]) - solve(Omega)[nodesInLayer, nodesInLayer]))
        }

      }
    }
    index_layer = sort(D2,index.return = T)$ix[1]
    sourceLayer = LayerToOrder[index_layer]
    LayerHOrder = c(LayerHOrder, sourceLayer)
    LayerToOrder = setdiff(c(1:layerNUM), LayerHOrder)

    index = index + 1
    layerOrder[[index]] = which(layerID == sourceLayer)
  }

  LayerHOrder = c(LayerHOrder,LayerToOrder)
  layerOrder[[layerNUM]] = which(layerID == LayerToOrder)

  return(list(layerOrder = layerOrder, LayerHOrder = LayerHOrder))
}
