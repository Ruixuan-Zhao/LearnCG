# LearnCG

## Introduction

The goal of this package is to study the Gaussian Chain Graph model. An efficient learning algorithm is built upon the identifiability conditions to fully recover the chain graph structure.

## Maintainer

Ruixuan Zhao (ruixuzhao2-c@my.cityu.edu.hk)

## Publication

Zhao, R., Zhang, H. & Wang, J. (2023+). Identifiability and Consistent Estimation of the Gaussian Chain Graph Model. Manuscript.

## Installation

Method 1: Run the following codes in R/Rstudio.

```
library(devtools)
devtools::install_github("Ruixuan-Zhao/LearnCG")
```

Method 2: Download the LearnCG_1.0.0.tar.gz and install from Package Archive File using RStudio.

## Usage

Detailed help functions are provided in the package for your reference.

### Toy Example

```
library(LearnCG)
## Generate a random chain graph with sample matrix
gen_data = simul.CG.Ex1(p.node = 50, n.sample = 500, p.un = 0.03,
                       p.hub = 0.2, p.hub.dir = 0.8, CG.seed = 1)
X = gen_data$X
## Outputs
est_res = Recover.CG(X = X, lambda1_n = 0.4*500^(-3/8),
                    lambda2_n = 0.8*500^(-3/8), kappa_n = 5.5*500^(-3/8),
                    nu_n = 3.8*500^(-1/4))
est_Omega = est_res$est_Omega
est_B = est_res$est_B
```
