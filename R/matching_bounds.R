cplex_MBs_problem <- function(Yt, Yc, m, kt=1, kc=1, constraints=list(), approximate=FALSE){

  nt = length(Yt)
  nc = length(Yc)

  tloc = gen_tloc(nt, nc)
  cloc = gen_cloc(nt, nc)

  cvec = as.numeric(sapply(Yt, "-", Yc))

  a1 = matrix(1, 1, nt*nc)
  b1 = m
  s1 = "E"

  a2 = matrix(0, nt, nt*nc)
  for (i in 1:nt){
    a2[i, tloc(i)] = 1
  }
  b2 = rep(kt, nt)
  s2 = rep("L", nt)

  a3 = matrix(0, nc, nt*nc)
  for (i in 1:nc){
    a3[i, cloc(i)] = 1
  }
  b3 = rep(kc, nc)
  s3 = rep("L", nc)

  if(approximate)
    vtype="C"
  else
    vtype="B"

  args = list(cvec=cvec,
              Amat=rbind(a1, a2, a3),
              bvec=c(b1, b2, b3),
              sense=c(s1, s2, s3),
              vtype=vtype,
              ub=1, lb=0)

  for(con in constraints){
    if(!is.null(con$Amat))
      args$Amat = rbind(args$Amat, con$Amat)
    if(!is.null(con$bvec))
      args$bvec = c(args$bvec, con$bvec)
    if(!is.null(con$sense))
      args$sense = c(args$sense, con$sense)
    if(!is.null(con$ub))
      args$ub = con$ub
  }
  args
}

matching_bounds <- function(Y, D, constraints, M, Kt, Kc, include_upper=TRUE,
                            include_lower=TRUE, approximate=TRUE,
                            approx_type="fast", solver="CPLEX", solver_control=list(),
                            return_solver_out=FALSE){

  nt <- sum(D==1)
  nc <- sum(D==0)
  message("Generating constraints...")
  generated_constraints = list()
  for (el in 1:length(constraints)){
    generated_constraints[[el]] = constraints[[el]](D)
  }

  message("Generating problem parameters...")
  args = cplex_MBs_problem(Y[D==1], Y[D==0], M, Kt, Kc, generated_constraints, approximate)

  if(solver=="CPLEX"){
    solver_res <- matching_bounds_cplex(args, nt, nc, Kt, Kc, include_upper,
                                        include_lower, approximate, approx_type,
                                        solver_control, return_solver_out)
  } else if(solver=="GLPK"){
    solver_res <- matching_bounds_glpk(args, nt, nc, Kt, Kc, include_upper,
                                        include_lower, approximate, approx_type,
                                        solver_control, return_solver_out)
  }else{
    stop("Solver must be one of 'GLPK' or 'CPLEX'")
  }

  res <- solver_res
  max_m <- solver_res$max_m
  min_m <- solver_res$min_m
  res$max_te_mb <- mean(Y[D==1][max_m[, 1]] - Y[D==0][max_m[, 2]])
  res$min_te_mb <- mean(Y[D==1][min_m[, 1]] - Y[D==0][min_m[, 2]])

  return(res)
}

matching_bounds_glpk <- function(args, nt, nc, Kt, Kc, include_upper=TRUE,
                                 include_lower=TRUE, approximate=TRUE,
                                 approx_type="fast", glpk_control=list(),
                                 return_solver_out=FALSE){

  dir <- sapply(args$sense, function(x){
    if(x == "L")
      return("<=")
    else if (x=="G")
      return(">=")
    else
      return("==")
  })
  args$sense <- dir

  if(length(args$lb)==1){
    lower <- rep(args$lb, length(args$cvec))
  }else{
    lower <- args$lb
  }

  if(length(args$ub)==1){
    upper <- rep(args$ub, length(args$cvec))
  }else{
    upper <- args$ub
  }

  args$bounds <- list(upper, lower)

  res <- list()

  if(include_upper){
    message("Computing Upper MB...")
    max_mb_l <- Rglpk::Rglpk_solve_LP(obj=args$cvec, mat=args$Amat,
                              rhs = args$bvec, dir = args$sense,
                              bounds = args$bounds,
                              types = args$vtype, max = TRUE,
                              control=glpk_control)
    if(!max_mb_l$status){
      warning("Solver failed to obtain a solution Upper MB")
    }

    message("Computing Integer approximation...")
    if(approximate){
      if(approx_type == "feasible"){
        max_mb <- approx_int_solution2_glpk(args, max_mb_l$solution, glpk_control)
        if(!max_mb$status){
          warning("Solver failed to obtain an integer solution for the Upper MB")
        }
      }else{
        max_mb <- approx_int_solution_glpk(nt, nc, Kt, Kc,
                                            max_mb_l$solution, glpk_control)
        if(!max_mb$status){
          warning("Solver failed to obtain an integer solution for the Upper MB")
        }
      }
    }

    max_m = which(matrix(max_mb$solution, nt, nc, byrow=T)==1, arr.ind=T)
    res$max_m <- max_m
  }

  if(include_lower){
    message("Computing Lower MB...")
    min_mb_l <- Rglpk::Rglpk_solve_LP(obj=args$cvec, mat=args$Amat,
                               rhs = args$bvec, dir = dir,
                               bounds = list(upper=upper, lower=lower),
                               types = args$vtype, max = FALSE,
                               control=glpk_control)
    if(!max_mb_l$status){
      warning("Solver failed to obtain a solution lower MB")
    }
    message("Computing Integer approximation...")
    if(approximate){
      if(approx_type == "feasible"){
        min_mb <- approx_int_solution2_glpk(args, min_mb_l$solution, glpk_control)
        if(!min_mb$status){
          warning("Solver failed to obtain an integer solution for the lower MB")
        }
      }else{
        min_mb <- approx_int_solution_glpk(nt, nc, Kt, Kc,
                                           min_mb_l$solution, glpk_control)
        if(!min_mb$status){
          warning("Solver failed to obtain an integer solution for the lower MB")
        }
      }
    }

    min_m = which(matrix(min_mb$solution, nt, nc, byrow=T)==1, arr.ind=T)
    res$min_m <- min_m
  }

  if(return_solver_out){
    res$solver_out_max <- max_mb_l
    res$solver_out_min <- min_mb_l
  }

  return(res)
}

matching_bounds_cplex <- function(args, nt, nc, Kt, Kc, include_upper=TRUE,
                                  include_lower=TRUE, approximate=TRUE,
                                  approx_type="fast", cplex_control=list(),
                                  return_cplex_out=FALSE){

  res <- list()

  if(include_upper){
    message("Computing Upper MB...")
    max_mb_l = do.call(Rcplex::Rcplex, args=c(args, objsense='max',
                                              control=list(cplex_control)))

    message("Computing Integer approximation...")
    if(approximate){
      if(approx_type == "feasible")
        max_mb <- approx_int_solution2(args, max_mb_l$xopt, cplex_control)
      else
        max_mb <- approx_int_solution(nt, nc, Kt, Kc,
                                      max_mb_l$xopt, cplex_control)
    }

    max_m = which(matrix(max_mb$xopt, nt, nc, byrow=T)==1, arr.ind=T)
    res$max_m <- max_m
  }

  if(include_lower){
    message("Computing Lower MB...")
    min_mb_l = do.call(Rcplex::Rcplex, args=c(args, objsense='min',
                                              control=list(cplex_control)))

    message("Computing Integer approximation...")
    if(approximate){
      if(approx_type == "feasible")
        min_mb <- approx_int_solution2(args, min_mb_l$xopt, cplex_control)
      else
        min_mb <- approx_int_solution(nt, nc, Kt, Kc,
                                        min_mb_l$xopt, cplex_control)
    }
    min_m = which(matrix(min_mb$xopt, nt, nc, byrow=T)==1, arr.ind=T)
    res$min_m <- min_m
  }

  if(return_cplex_out){
    res$cplex_out_max <- max_mb_l
    res$cplex_out_min <- min_mb_l
  }

  return(res)
}
