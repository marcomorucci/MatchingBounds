gen_tloc <- function(nt, nc){
  function(i){
    (nc*(i-1)+1):(nc*i)
  }
}

gen_cloc <- function(nt, nc){
  function(i){
    (0:(nt-1))*nc + i
  }
}

gen_ploc <- function(nt, nc){
  tloc <- gen_tloc(nt, nc)
  cloc <- gen_cloc(nt, nc)
  ploc <- function(i, j){
    intersect(tloc(i), cloc(j))
  }
}

approx_int_solution <- function(nt, nc, kt, kc, t_strict, c_strict, coef, cplex_control=list()){
  tloc <- gen_tloc(nt, nc)
  cloc <- gen_cloc(nt, nc)

  a2 <- matrix(0, nt, nt*nc)
  for (i in 1:nt){
    a2[i, tloc(i)] <- 1
  }
  b2 <- rep(kt, nt)
  if (!t_strict)
    s2 <- rep("L", nt)
  else
    s2 <- rep("E", nt)

  a3 <- matrix(0, nc, nt*nc)
  for (i in 1:nc){
    a3[i, cloc(i)] <- 1
  }
  b3 <- rep(kc, nc)
  if (!c_strict)
    s3 <- rep("L", nc)
  else
    s3 <- rep("E", nc)

  Amat <- rbind(a2, a3)
  bvec <- c( b2, b3)
  sense <- c(s2, s3)
  cvec <- coef
  ub <- Inf
  vtype <- rep("B", nt*nc)

  Rcplex::Rcplex(cvec, Amat, bvec, sense=sense, ub=ub,
                    objsense='max', vtype = vtype, control=cplex_control)
}


approx_int_solution_glpk <- function(nt, nc, kt, kc, t_strict, c_strict, coef, glpk_control=list()){
  tloc <- gen_tloc(nt, nc)
  cloc <- gen_cloc(nt, nc)

  a2 <- matrix(0, nt, nt*nc)
  for (i in 1:nt){
    a2[i, tloc(i)] <- 1
  }
  b2 <- rep(kt, nt)
  if (!t_strict)
    s2 <- rep("<=", nt)
  else
    s2 <- rep("=", nt)

  a3 <- matrix(0, nc, nt*nc)
  for (i in 1:nc){
    a3[i, cloc(i)] <- 1
  }
  b3 <- rep(kc, nc)
  if (!c_strict)
    s3 <- rep("<=", nc)
  else
    s3 <- rep("=", nc)

  Amat <- rbind(a2, a3)
  bvec <- c( b2, b3)
  sense <- c(s2, s3)
  vtype <- rep("B", nt*nc)

  Rglpk::Rglpk_solve_LP(obj=coef, mat=args$Amat,
                 rhs = bvec, dir = sense,
                 types = vtype, max = TRUE,
                 control=glpk_control)
}

approx_int_solution2 <- function(args, coef, cplex_control=list()){
  Rcplex::Rcplex(Amat=args$Amat,
                        bvec=args$bvec, cvec=coef, ub=args$ub,
                        vtype = "B", n=1, sense=args$sense,
                        objsense='max', control = cplex_control)
}

approx_int_solution2_glpk <- function(args, coef, glpk_control=list()){
  Rglpk::Rglpk_solve_LP(obj=coef, mat=args$Amat,
                 rhs = args$bvec, dir = args$sense,
                 bounds = args$bounds,
                 types = rep("B", length(args$cvec)), max = TRUE,
                 control=glpk_control)
}


flip_constraints <- function(nt, nc, constraints){
  cloc<-gen_cloc(nt, nc)
  flipped_constraints <- list()
  for (con in names(constraints)){
    flipped_constraints[[con]] <- constraints[[con]]
    if(!is.null(constraints[['con']]$Amat))
      flipped_constraints[[con]] <- constraints[[con]]$Amat[, c(t(sapply(1:nc, cloc)))]
    if(!is.null(constraints[[con]]$ub))
      flipped_constraints[[con]]$ub <- constraints[[con]]$ub[c(t(sapply(1:nc, cloc)))]
  }
  flipped_constraints
}
