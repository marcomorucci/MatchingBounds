cplex_MBs_problem = function(Yt, Yc, m, kt=1, kc=1,
                                constraints=list(), approximate=FALSE){
  
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


matching_bounds <- function(Y, D, constraints, M, Kt, Kc, approximate=TRUE){
    generated_constraints = list()
    for (el in 1:length(constraints)){
        generated_constraints[[el]] = constraints[[el]](D)
    }
    args = cplex_Z_problem(Y[D==1], Y[D==0], M, Kt, Kc, generated_constraints, approximate)

    max_mb = do.call(Rcplex::Rcplex, args=c(args, objsense='max', control=list(list(trace=0))))
    if(approximate)
      max_mb = approx_int_solution2(args, max_mb$xopt, cplex_control)
    max_m = which(matrix(max_mb$xopt, nt, nc, byrow=T)==1, arr.ind=T)
    max_te_mb = mean(Yt[max_m[, 1]] - Yc[max_m[, 2]])
    
    min_mb = do.call(Rcplex::Rcplex, args=c(args, objsense='min', control=list(list(trace=0))))
    min_mb = approx_int_solution2(args, min_mb$xopt, cplex_control)
    if(approximate)
      max_mb = approx_int_solution2(args, max_mb$xopt, cplex_control)
    min_m = which(matrix(min_mb$xopt, nt, nc, byrow=T)==1, arr.ind=T)
    min_te_mb = mean(Yt[min_m[, 1]] - Yc[min_m[, 2]])
    
    return(list(max_m, min_m, max_te_mb, min_te_mb))
}
