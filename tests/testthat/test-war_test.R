library(dplyr)

dta <- read.csv("tests/adi.csv")

Y <- dta$idp
D <- dta$actualviolence
X <- select(dta, industry, income,
            crop_cattle_seized, landseized_self, homedestr_cause,
            dstryd_indstry, communityorg, motorroad, maoists,
            totalchildren, education, newgender, age) %>%
  mutate(income=scale(income), totalchildren=scale(totalchildren),
         education=scale(education), age=scale(age)) %>%
  as.matrix(.)

psc <- predict(glm(D ~ X, family=binomial(link='logit')), type='response')
dmat = abs(t(sapply(psc, "-", psc)))
l2_dmat <- apply(X, 1, function(x)  sqrt(colSums((t(X) - x)^2)))

# ATC
Yt <- Y[D==1]
Yc <- Y[D==0]
Xt <- X[D==1, ]
Xc <- X[D==0, ]
nt <- sum(D)
nc <- sum(1-D)
M <- nc
Kt <- nc
Kc <- 1

#1nn psc match
Ym <- rep(NA, nc)
dm <- rep(NA, nc)
Xm <- matrix(NA, nc, ncol(X))
pscm <- rep(NA, nc)
for(j in 1:nc){
  dist <- dmat[D==0, D==1][j, ] #sqrt(rowSums((Xt - Xc[j, ])^2))
  mj <- order(dist)[1]
  Ym[j] <- Yt[mj]
  Xm[j, ] <- Xt[mj, ]
  dm[j] <- dmat[D==0, D==1][j, mj]
  pscm[j] <- psc[D==0][mj]
}

nn_t <- t.test(Ym, Yc)
Xnn <- rbind(Xm, Xc)
Dnn <- rep(c(1,0), each=M)
Ynn <- c(Ym, Yc)
lmnn <- lm(Ynn ~ Dnn + Xnn)


approximate = TRUE
cplex_control = list(trace = 1, round = 1, tilim = 60*60)
mean_tol <- abs(colMeans(Xm - Xc))
var_tol <- abs(colMeans(Xm^2 - Xc^2))
skw_tol <- abs(colMeans(Xm^3 - Xc^3))
dm_mat <- matrix(rep(dm, each=nt), nt, nc, byrow=FALSE)
epsilons <- seq(1, 2, by=0.1)
res <- NULL
for(eps in epsilons){
  constraints = list(caliper=constraint_caliper(dmat, dm_mat * eps ),
                     mean=constraint_moment(X, M, 1, mean_tol * eps + 1e-5),
                     var=constraint_moment(X, M, 2, var_tol * eps + 1e-5),
                     skw=constraint_moment(X, M, 3, skw_tol * eps + 1e-5))

  generated_constraints = list()
  for (el in 1:length(constraints)){
    generated_constraints[[el]] = constraints[[el]](D)
  }

  args = cplex_MBs_problem(Yt, Yc, M, Kt, Kc, generated_constraints, approximate)
  max_mb = do.call(Rcplex::Rcplex, args=c(args, objsense='max', control=list(cplex_control)))
  max_mb = approx_int_solution2(args, max_mb$xopt, list(cplex_control))
  max_m = which(matrix(max_mb$xopt, nt, nc, byrow=T)==1, arr.ind=T)

  min_mb = do.call(Rcplex::Rcplex, args=c(args, objsense='min', control=list(cplex_control)))
  min_mb = approx_int_solution2(args, min_mb$xopt, list(cplex_control))
  min_m = which(matrix(min_mb$xopt, nt, nc, byrow=T)==1, arr.ind=T)

  res <- rbind(res, c(eps=eps, bd=1, max_m[, 1]), c(eps=eps, bd=0, min_m[, 1]))
}
