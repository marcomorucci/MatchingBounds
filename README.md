
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MatchingBounds

<!-- badges: start -->
<!-- badges: end -->

Matching Bounds implements a set of integer-programming based
diagnostics aimed at exploring whether multiple matched sets exists in
the same dataset that produce different estimates of the causal
quantities of interest.

## Installation

You can install the released version of MatchingBounds from
[Github](https://github.com/marcomorucci/MatchingBounds) with:

``` r
devtools::install_github("https://github.com/marcomorucci/MatchingBounds")
```

## Example

This is a simple example of applying matching bounds to the subset of
the [Lalonde data](https://users.nber.org/~rdehejia/data/.nswdata2.html)
available in the MatchIt library. This data consists of a set of treated
and control units to which a treatment was administered experimentally
(i.e., by choosing recipients randomly) as well as a set of control
units to which the treatment was not administered at random
(self-selection). The aim is to estimate the ATT of the treatment by
matching the treated experimental units to the non-experimental
controls.

First, we load the data:

``` r
library(MatchIt)
data(lalonde)

lalonde$black = lalonde$race=="black"
lalonde$hisp = lalonde$race=="hispan"

cvars <- c("age", "educ", "black", "hisp", "married", "nodegree", "re74", "re75")
```

Now we set up some hyperparameters of our matching methods. We set the
total number of matches to the number of treated units and we disallow
re-use of treated and control units: this will let us estimate the ATT
by one-to-one matching each treated unit to one control unit with no
replacement.

``` r
# This is the total number of matches we will make.
M <- sum(lalonde$treat)
# This is how many times we will allow re-use of treated units
Kt <- 1
# This is how many times we will allow re-use of control units
Kc <- 1
```

Then we create a constraint on the match quality in that the marginal
distributions of the covariates in the matched sets found by MBs must
have the same difference in means as those in the experimental sets of
the Lalonde data.

``` r
# These are the differences in covariate means of the treated experimental sample - the control experimental sample.
tols <- c(age=0.76237006, educ=0.25748441,   black=0.01632017,
          hisp=0.04823285,  married=0.03534304, nodegr=0.12650728,
          re74=11.45281538, re75=265.14638896)


constraints = list(mean=constraint_moment(lalonde[, cvars], lalonde$treat, M, 1, tols))
```

Finally, we run the matching bounds algorithm.

``` r
# Solve using CPLEX
res <- matching_bounds(lalonde$re78, lalonde$treat, constraints,
                       M, Kt, Kc, approximate=TRUE,
                       include_upper = TRUE,
                       include_lower = TRUE, approx_type = 'fast', solver = 'CPLEX',
                       solver_control = list(trace = 0, round = 1, tilim = 60*60),
                       return_solver_out = TRUE)
# If using GLPK 
# res <- matching_bounds(lalonde$re78, lalonde$treat, constraints,
#                        M, Kt, Kc, approximate=TRUE,
#                        include_upper = TRUE,
#                        include_lower = TRUE, approx_type = 'fast', solver = 'GLPK',
#                        solver_control = list(verbose=TRUE, presolve=TRUE, tm_limit= 60*60*1000),
#                        return_solver_out = TRUE)
```

The largest ATT that can be found among matches that satisfy the balance
constraint is 504.13 and the smallest is 73.32.

We can also run MBs many times with progressively looser constraints
using the progressive functions in the package. Here we run MBs six
times total, each time relaxing the constraints on the means by 10% of
its initial value.

``` r

means <- progressive_constraint_moment(lalonde[, cvars], lalonde$treat, M, 1, tols,
                                       rf_step=0.1, end_rf=0.6)

progressive_constraints <- list(means=means)

prog_res <- progressive_matching_bounds(lalonde$re78, lalonde$treat, progressive_constraints,
                                   M, Kt, Kc, approximate=TRUE,
                                   include_upper = TRUE,
                                   include_lower = TRUE, approx_type = 'fast', solver = 'CPLEX',
                                   solver_control = list(trace = 0, round = 1, tilim = 60*60),
                                   return_solver_out = TRUE)
```

Results can be viewed as follows:

``` r
summary(prog_res)
#>      TE: UB     TE: LB N matched T: UB N matched C: UB N matched T: LB
#> 1  504.1336   73.31937             126             126             121
#> 2 1179.4261 1063.29831             124             124             121
#> 3 1048.6568  816.14680             123             123             122
#> 4 2150.7023  735.70870             121             121             121
#> 5  821.4618 -477.00164             123             123             125
#> 6 1105.3041   17.96618             122             122             123
#>   N matched C: LB means RF
#> 1             121      0.0
#> 2             121      0.1
#> 3             122      0.2
#> 4             121      0.3
#> 5             125      0.4
#> 6             123      0.5
```
