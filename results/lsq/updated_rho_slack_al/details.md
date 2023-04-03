Ran ALBO without L-BFGS-B optimising the EI acquisition function at each iteration on
the toy problem, with the objective treated as being *unknown* and modelled using a GP.

NOTE - This run used the correct method for initialising the penalty parameter.

Code for run can be found below.

```
library(laGP)
library(jsonlite)
library(glue)
library(R.utils)

new_rho_update <- function (obj, C, equal, init = 10, ethresh = 0.1)
{
    print("Thomas Code")
    if (init > length(obj))
        init <- length(obj)
    if (length(obj) != nrow(C))
        stop("length(obj) != nrow(C)")
    v <- apply(C, 1, function(x) {
        all(x[!equal] <= 0) && all(abs(x[equal]) < ethresh)
    })
    Civ <- C[!v, , drop = FALSE]
    if (nrow(Civ) == 0)
      return (1/2)
    valid <- Civ <= 0
    inequality <- matrix(!equal, nrow=nrow(Civ), ncol=length(!equal), byrow=TRUE)
    Civ[valid & inequality] <- 0  # Only replace *inequality constraints* <= 0 with 0
    Civs2 <- rowSums(Civ^2)
    cm2 <- min(Civs2)
    f <- obj[1:init]
    fv <- f[v[1:init]]
    if (length(fv) > 0)
        fm <- min(fv)
    else fm <- median(f)
    rho <- cm2/(2 * abs(fm))
    return(rho)
}

assignInNamespace("auto.rho", new_rho_update, ns = "laGP")

## toy function returning linear objective evaluations and
## non-linear constraints
aimprob <- function(X, known.only = FALSE)
{
  if(is.null(nrow(X))) X <- matrix(X, nrow=1)
  f <- rowSums(X)
  if(known.only) return(list(obj=f))
  c1 <- 1.5-X[,1]-2*X[,2]-0.5*sin(2*pi*(X[,1]^2-2*X[,2]))
  c2 <- rowSums(X^2)-1.5
  return(list(obj=f, c=cbind(c1,c2)))
}

## set bounding rectangle for adaptive sampling
B <- matrix(c(rep(0,2),rep(1,2)),ncol=2)

for(x in 1:100) {
  ## run ALBO
  out <- optim.auglag(aimprob, B, start=5, end=50, slack=TRUE, fhat=TRUE, lambda=0)
  write_json(out, glue("results/lsq/updated_rho_slack_al/data/run_{x}_results.json"), digits=NA)
}



```


