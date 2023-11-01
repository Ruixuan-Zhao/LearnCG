################################################################################
######################## Generate Chain Graph ##################################
################################################################################

#' Function to generate random chain graph (Example 1)
#'
#' @description The function to generate random chain graph.
#'
#' @details Particularly, we randomly connect each pair of nodes by undirected
#' edge with a certain probability
#' and read off multiple chain components. Then, a causal ordering of the chain
#' components is randomly assigned. For each chain component, we randomly select
#' the nodes as hubs with a certain probability, and let each hub node points to
#' the nodes in lower chain components with a certain probability.
#'
#'
#'
#' @usage simul.CG.Ex1(p.node, n.sample, p.un, p.hub, p.hub.dir, CG.seed)
#'
#' @param p.node Number of nodes
#' @param n.sample Sample size
#' @param p.un The probability of connecting two nodes with undirected edge
#' @param p.hub The probability of the number of hub nodes in each chain component
#' @param p.hub.dir The probability of directed edges pointing from a hub node to the nodes in the lower chain components
#' @param CG.seed The random number for generating random chain graph
#'
#' @return A list including
#' \describe{
#' \item{Omega}{The precision matrix of noise}
#' \item{B}{The coefficient matrix corresponding to the directed edges}
#' \item{X}{The data matrix}
#' }
#'
#' @export
#'
#' @import igraph mvtnorm
#'
#' @importFrom stats rbinom runif
#'
#' @examples
#' simul.CG.Ex1(p.node = 10, n.sample = 200, p.un = 0.1, p.hub = 0.5, p.hub.dir = 0.5, CG.seed = 1)
#'
#'

################################################################################
############################### Example 1 ######################################
################################################################################
# Remark: This is corresponding to the Example 2 in paper.

simul.CG.Ex1 = function(p.node, n.sample, p.un, p.hub, p.hub.dir, CG.seed){
  ## ---------------------------------------------------------------------------
  ## The name of the function: simul.CG.Ex1
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to generate random chain graph.
  ## ---------------------------------------------------------------------------
  ## Input:
  ## @ p.node: number of nodes
  ## @ n.sample: sample size
  ## @ p.un: the probability of connecting two nodes with undirected edge
  ## @ p.hub: the probability of the number of hub nodes in each chain component
  ## @ p.hub.dir: the probability of directed edges pointing from a hub node to
  ##            the nodes in the lower chain components
  ## @ CG.seed: the random number for generating random chain graph
  ## ---------------------------------------------------------------------------
  ## Output:
  ## @ Omega: precision matrix of noise corresponding to the undirected edges
  ## @ B: coefficient matrix corresponding to the directed edges
  ## @ X: sample matrix
  ## ---------------------------------------------------------------------------

  set.seed(CG.seed)

  Omega = matrix(0,p.node,p.node)
  B = matrix(0,p.node,p.node)

  ### Generate Chain graph
  # generate undirected edges
  Omega.pre = Omega
  w = which(lower.tri(Omega.pre))
  Omega.pre[w] = rbinom(n=length(w), size = 1, prob = p.un)*sample(c(-1,1),length(w),replace=T)*runif(n=length(w),min=0.5,max=1.5)
  Omega.pre = Omega.pre + t(Omega.pre)
  # permute to diagonal block
  unG = igraph::graph_from_adjacency_matrix(Omega.pre,mode = "undirected", weighted = TRUE)
  ccINF = igraph::components(unG)
  ccID = sort(ccINF$membership, index.return = TRUE)$ix

  Omega = Omega.pre[ccID,ccID]

  m.cc = ccINF$no

  # generate chain component
  randomCC = cbind(1:p.node, sort(ccINF$membership))
  cc.list = sapply(1:m.cc, function(x) randomCC[randomCC[,2]==x,1],simplify = F)

  diag(Omega) = colSums(abs(Omega)) + 0.1

  # directed edges across chain components
  if (m.cc > 1){
    lower.cc = 1:p.node
    for (i in 1:(m.cc-1)){
      current.cc = cc.list[[i]]
      lower.cc = setdiff(lower.cc, cc.list[[i]])
      set.seed(i + CG.seed)
      hub.nodes = rbinom(n = length(current.cc), size = 1, prob = p.hub)
      if (length(hub.nodes) > 0){
        index.hub.nodes = current.cc[hub.nodes==1]
        B[lower.cc, index.hub.nodes] = matrix(rbinom(n = length(index.hub.nodes)*length(lower.cc), size = 1, prob = p.hub.dir)*sample(c(-1,1), length(index.hub.nodes)*length(lower.cc), replace = T)*runif(n = length(index.hub.nodes)*length(lower.cc), min = 0.5, max = 1.5), length(lower.cc), length(index.hub.nodes))
      }
    }
  }


  ### Generate data
  Ip = diag(p.node)
  Theta = t(Ip - B)%*%Omega%*%(Ip - B)
  Sigma = solve(Theta)
  X = mvtnorm::rmvnorm(n.sample, mean = rep(0, p.node), sigma = Sigma)

  return(list(Omega = Omega, B = B, X = X))
}





#' Function to generate random chain graph (Example 2)
#'
#' @description The function to generate random chain graph.
#'
#' @details  This random chain graph is also called two-layer Gaussian graphical model.
#' Particularly, we randomly assign all the nodes into two layers. Within each
#' layer, we randomly connect each pair of nodes by an undirected edge with a
#' certain probability. Then, we generate the directed edges from
#' nodes in one layer (upper layer) to nodes in another layer (lower layer) with a certain
#'  probability.
#'
#'
#' @usage simul.Two.Layer.CG(p.node, n.sample, p.upper, p.un, p.dir, CG.seed)
#'
#' @param p.node Number of nodes
#' @param n.sample Sample size
#' @param p.upper The probability of assigning the nodes to the upper layer
#' @param p.un The probability of connecting two nodes with undirected edge in two layers
#' @param p.dir The probability of directed edges pointing from the node in the upper layer to the node in the lower layer
#' @param CG.seed The random number for generating random chain graph
#'
#' @return A list including
#' \describe{
#' \item{Omega}{The precision matrix of noise}
#' \item{B}{The coefficient matrix corresponding to the directed edges}
#' \item{X}{The data matrix}
#' }
#'
#' @export
#'
#' @import igraph mvtnorm
#'
#' @importFrom stats rbinom runif
#'
#' @examples
#' \donttest{
#' simul.Two.Layer.CG(p.node=20, n.sample=100, p.upper=0.4, p.un=0.05, p.dir=0.5, CG.seed=2)
#' }
#'




################################################################################
############################### Example 2 ######################################
################################################################################
# Remark: This is corresponding to the Example 1 in paper.

### two-layer gaussian chain graph
simul.Two.Layer.CG = function(p.node, n.sample, p.upper, p.un, p.dir, CG.seed){
  ## ---------------------------------------------------------------------------
  ## The name of the function: simul.Two.Layer.CG
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to generate random chain graph.
  ## ---------------------------------------------------------------------------
  ## Input:
  ## @ p.node: number of nodes
  ## @ n.sample: sample size
  ## @ p.upper: the probability of assigning the nodes to the upper layer
  ## @ p.un: the probability of connecting two nodes with undirected edge in two layers
  ## @ p.dir: the probability of directed edges pointing from the node in the upper layer to
  ##            the node in the lower layer
  ## @ CG.seed: the random number for generating random chain graph
  ## ---------------------------------------------------------------------------
  ## Output:
  ## @ Omega: precision matrix of noise corresponding to the undirected edges
  ## @ B: coefficient matrix corresponding to the directed edges
  ## @ X: sample matrix
  ## ---------------------------------------------------------------------------

  set.seed(CG.seed)

  Omega = matrix(0,p.node,p.node)
  B = matrix(0,p.node,p.node)

  upper.cc.cut = ceiling(p.node*p.upper)

  ### Generate chain graph
  ## Generate undirected edges
  Omega.upper.pre = matrix(0,upper.cc.cut, upper.cc.cut)
  w.upper = which(lower.tri(Omega.upper.pre))
  Omega.upper.pre[w.upper] = rbinom(n=length(w.upper), size = 1, prob = p.un)*sample(c(-1,1),length(w.upper),replace=T)*runif(n=length(w.upper),min=0.5,max=1.5)
  Omega.upper = Omega.upper.pre + t(Omega.upper.pre)

  Omega.lower.pre = matrix(0,p.node - upper.cc.cut, p.node - upper.cc.cut)
  w.lower = which(lower.tri(Omega.lower.pre))
  Omega.lower.pre[w.lower] = rbinom(n=length(w.lower), size = 1, prob = p.un)*sample(c(-1,1),length(w.lower),replace=T)*runif(n=length(w.lower),min=0.5,max=1.5)
  Omega.lower = Omega.lower.pre + t(Omega.lower.pre)

  Omega[1:upper.cc.cut, 1:upper.cc.cut] = Omega.upper
  Omega[(upper.cc.cut + 1):p.node, (upper.cc.cut + 1):p.node] = Omega.lower
  diag(Omega) = colSums(abs(Omega)) + 0.1

  ## Generate directed edges
  B[(upper.cc.cut + 1):p.node, 1:upper.cc.cut] = matrix(rbinom(n=(p.node - upper.cc.cut)*upper.cc.cut, size = 1, prob = p.dir)*sample(c(-1,1),(p.node - upper.cc.cut)*upper.cc.cut,replace=T)*runif(n=(p.node - upper.cc.cut)*upper.cc.cut,min=0.5,max=1.5), p.node - upper.cc.cut, upper.cc.cut)

  ### Generate sample matrix
  ### Generate data
  Ip = diag(p.node)
  Theta = t(Ip-B)%*%Omega%*%(Ip-B)
  Sigma = solve(Theta)
  X = mvtnorm::rmvnorm(n.sample,mean=rep(0,p.node),sigma=Sigma)


  return(list(Omega = Omega, B = B, X = X))
}




