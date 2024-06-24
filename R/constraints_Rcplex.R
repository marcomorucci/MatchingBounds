#' Define a distance caliper constraint
#'
#' @param dist_matrix a N x N matrix of distances between units.
#' @param caliper either: a Nt x Nc matrix of maximum allowed distances between treated and control units or a scalar denoting the max allowed distance between any pair of units.
#'
#' @return A function to be used in conjunction with one of the MB computations
#' @export
#'
#' @examples
constraint_caliper <- function(dist_matrix, caliper){

  if(is.null(dim(caliper))){
    add_constraint_caliper = function(D){
      dmat = dist_matrix[D==1, D==0]
      ub = as.numeric(t(dmat) <= caliper)
      list(ub=ub)
    }
  }else{
      add_constraint_caliper = function(D){
        dmat = dist_matrix[D==1, D==0]
        ub = as.numeric(t(dmat) <= t(caliper))
        list(ub=ub)
    }
  }
  add_constraint_caliper
}


#' Defines a moment distance constraint
#'
#' @param X matrix of matching covariates.
#' @param m Integer. total number of matches to be made.
#' @param mom Integer. moment to constrain.
#' @param eps Either a vector of numeric tolerances with length(eps) == ncol(X) or a scalar defining a global tolerance for all covariates.
#'
#' @return A function to be used in conjunction with one of the MB computations
#' @export
#'
#' @examples
constraint_moment <- function(X, m, mom, eps){
  add_constraint_moment = function(D){
    Xt = X[D==1, ]
    Xc = X[D==0, ]

    nt = nrow(Xt)
    nc = nrow(Xc)
    p = ncol(Xt)

    if(length(eps) == 1)
      eps = rep(eps, p)

    a = matrix(0, p, nt*nc)
    for (j in 1:p){
      a[j, ] = as.numeric(sapply(Xt[, j]**mom, "-", Xc[, j]**mom))/m
    }
    b = eps
    s = rep("L", p)

    Amat = rbind(a, -a)
    bvec = c(b, b)
    sense = c(s, s)
    list(Amat=Amat, bvec=bvec, sense=sense)
  }
}

#' Defines a quantile distance constraint
#'
#' @param X matrix of matching covariates.
#' @param m Integer. total number of matches to be made.
#' @param quantiles List. A list of vectors with length(quantiles) == ncol(X). Each vector should specify the values of the Nq quantiles of each covariate that the user wish to constrain.
#' @param eps Either a vector of numeric tolerances with as many elements as the total elements in quantiles or a scalar defining a global tolerance for all quantiles of all covariates.
#'
#' @return A function to be used in conjunction with one of the MB computations
#' @export
#'
#' @examples
constraint_quantile <- function(X, m, quantiles, eps){

  add_constraint_quantile = function(D){
    Xt = X[D==1, ]
    Xc = X[D==0, ]
    nt = nrow(Xt)
    nc = nrow(Xc)
    p = ncol(Xt)

    a = matrix(0, sum(sapply(quantiles, length)), nt*nc)
    i = 1
    for (j in 1:p){
      for(q in quantiles[[j]]){
        a[i, ] = as.numeric(sapply((Xt[, j] <= q), "-", (Xc[, j] <= q)))/m
        i = i + 1
      }
    }
    if (length(eps)==1)
      b = rep(eps, sum(sapply(quantiles, length)))
    else
      b = unlist(eps)
    s = rep("L", sum(sapply(quantiles, length)))


    Amat = rbind(a, -a)
    bvec = c(b, b)
    sense = c(s, s)
    list(Amat=Amat, bvec=bvec, sense=sense)
  }
}
