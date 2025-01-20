#' Define a distance caliper constraint
#'
#' @param dist_matrix EITEHR a N x N matrix of distances between units. Or a Nt x Nc matrix of distances between treated and control units. If the former the the argument D must also be provided.
#' @param caliper either: a Nt x Nc matrix of maximum allowed distances between treated and control units or a scalar denoting the max allowed distance between any pair of units.
#' @param D binary treatment indicator vector if not left NULL, then dist_matrix is assumed to be NxN
#'
#' @return A function to be used in conjunction with one of the MB computations
#' @export
#'
constraint_caliper <- function(dist_matrix, caliper, D=NULL){

  if(is.null(D))
    dmat <- dist_matrix
  else
    dmat = dist_matrix[D==1, D==0]

  if(is.null(dim(caliper))){
    ub = as.numeric(t(dmat) <= caliper)
  }else{
    ub = as.numeric(t(dmat) <= t(caliper))
  }
  list(ub=ub)
}


#' Defines a progressively-relaxed caliper distance constraint
#'
#' @param dist_matrix EITEHR a N x N matrix of distances between units. Or a Nt x Nc matrix of distances between treated and control units. If the former the the argument D must also be provided.
#' @param caliper either: a Nt x Nc matrix of maximum allowed distances between treated and control units or a scalar denoting the max allowed distance between any pair of units.
#' @param D binary treatment indicator vector if not left NULL, then dist_matrix is assumed to be NxN
#' @param init_rf The initial relaxation factor (a proportion > 0) of the constraints
#' @param end_rf The final relaxation factor that should be reached from the initial in nsteps steps. If provided then rf_step or nsteps can be automatically computed depending on which of the two is provided.
#' @param rf_step The amount to increase the relaxation factor by each step. If provided then either end_rf or nsteps can be automatically computed depending on which of the two is provided.
#' @param nsteps The total number of times MBs should be estimated. f provided then either end_rf or rf_step can be automatically computed depending on which of the two is provided.
#'
#' @return A function to be used in conjunction with one of the MB computations
#' @export
#'
progressive_constraint_caliper <- function(dist_matrix, caliper, D=NULL, init_rf=0.0, end_rf=NULL,
                                          rf_step=NULL, nsteps=NULL){

  if (is.null(end_rf) + is.null(rf_step) + is.null(nsteps) > 1){
    stop("At least two of end_rf, rf_step, nsteps must be specified")
  }

  if (is.null(rf_step)){
    rf_step <- (end_rf - init_rf)/nsteps
  }
  if(is.null(nsteps)){
    nsteps <- (end_rf - init_rf)/rf_step
  }
  cons <- list()
  rf <- init_rf
  for(step in 1:nsteps){
    cons[[step]] <- c(constraint_caliper(dist_matrix, caliper + caliper * rf, D),
                      list(rf=rf))
    rf <- rf + rf_step
  }
  cons
}


#' Defines a moment distance constraint
#'
#' @param X matrix of matching covariates.
#' @param D Binary treatment indicator vector
#' @param m Integer. total number of matches to be made.
#' @param mom Integer. moment to constrain.
#' @param eps Either a vector of numeric tolerances with length(eps) == ncol(X) or a scalar defining a global tolerance for all covariates.
#'
#' @return A function to be used in conjunction with one of the MB computations
#' @export
#'
constraint_moment <- function(X, D, m, mom, eps){
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

#' Defines a progressively-relaxed moment distance constraint
#'
#' @param X matrix of matching covariates.
#' @param D Binary treatment indicator vector
#' @param m Integer. total number of matches to be made.
#' @param mom Integer. moment to constrain.
#' @param eps Either a vector of numeric tolerances with length(eps) == ncol(X) or a scalar defining a global tolerance for all covariates.
#' @param init_rf The initial relaxation factor (a proportion > 0) of the constraints
#' @param end_rf The final relaxation factor that should be reached from the initial in nsteps steps. If provided then rf_step or nsteps can be automatically computed depending on which of the two is provided.
#' @param rf_step The amount to increase the relaxation factor by each step. If provided then either end_rf or nsteps can be automatically computed depending on which of the two is provided.
#' @param nsteps The total number of times MBs should be estimated. f provided then either end_rf or rf_step can be automatically computed depending on which of the two is provided.
#'
#' @return A function to be used in conjunction with one of the MB computations
#' @export
#'
progressive_constraint_moment <- function(X, D, m, mom, eps, init_rf=0.0, end_rf=NULL,
                                          rf_step=NULL, nsteps=NULL){

  if (is.null(end_rf) + is.null(rf_step) + is.null(nsteps) > 1){
    stop("At least two of end_rf, rf_step, nsteps must be specified")
  }

  if (is.null(rf_step)){
    rf_step <- (end_rf - init_rf)/nsteps
  }
  if(is.null(nsteps)){
    nsteps <- (end_rf - init_rf)/rf_step
  }
  cons <- list()
  rf <- init_rf
  for(step in 1:nsteps){
    cons[[step]] <- c(constraint_moment(X, D, m, mom, eps + eps * rf),
                      list(rf=rf))
    rf <- rf + rf_step
  }
  cons
}

#' Defines a quantile distance constraint
#'
#' @param X matrix of matching covariates.
#' @param D Binary treatment indicator vector
#' @param m Integer. total number of matches to be made.
#' @param quantiles List. A list of vectors with length(quantiles) == ncol(X). Each vector should specify the values of the Nq quantiles of each covariate that the user wish to constrain.
#' @param eps Either a vector of numeric tolerances with as many elements as the total elements in quantiles or a scalar defining a global tolerance for all quantiles of all covariates.
#'
#' @return A function to be used in conjunction with one of the MB computations
#' @export
#'
constraint_quantile <- function(X, D, m, quantiles, eps){
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




#' Defines a progressively-relaxed quantile distance constraint
#'
#' @param X matrix of matching covariates.
#' @param D Binary treatment indicator vector
#' @param m Integer. total number of matches to be made.
#' @param quantiles List. A list of vectors with length(quantiles) == ncol(X). Each vector should specify the values of the Nq quantiles of each covariate that the user wish to constrain.
#' @param eps Either a vector of numeric tolerances with as many elements as the total elements in quantiles or a scalar defining a global tolerance for all quantiles of all covariates.
#'
#' @param init_rf The initial relaxation factor (a proportion > 0) of the constraints
#' @param end_rf The final relaxation factor that should be reached from the initial in nsteps steps. If provided then rf_step or nsteps can be automatically computed depending on which of the two is provided.
#' @param rf_step The amount to increase the relaxation factor by each step. If provided then either end_rf or nsteps can be automatically computed depending on which of the two is provided.
#' @param nsteps The total number of times MBs should be estimated. f provided then either end_rf or rf_step can be automatically computed depending on which of the two is provided.
#'
#' @return A function to be used in conjunction with one of the MB computations
#' @export
#'
progressive_constraint_quantile <- function(X, D, m, quantiles, eps, init_rf=0.0, end_rf=NULL,
                                          rf_step=NULL, nsteps=NULL){

  if (is.null(end_rf) + is.null(rf_step) + is.null(nsteps) > 1){
    stop("At least two of end_rf, rf_step, nsteps must be specified")
  }

  if (is.null(rf_step)){
    rf_step <- (end_rf - init_rf)/nsteps
  }
  if(is.null(nsteps)){
    nsteps <- (end_rf - init_rf)/rf_step
  }
  cons <- list()
  rf <- init_rf
  for(step in 1:nsteps){
    cons[[step]] <- c(constraint_quantile(X, D, m, quantiles, eps + eps * rf),
                      list(rf=rf))
    rf <- rf + rf_step
  }
  cons
}

