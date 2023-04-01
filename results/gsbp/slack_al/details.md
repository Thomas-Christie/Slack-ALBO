Ran ALBO without L-BFGS-B optimising the EI acquisition function at each iteration on
the GSBP, with the objective treated as being *unknown* and modelled using a GP.

NOTE - This run used the author's (faulty) code for initialising penalty term.

Code for run can be found below.

```
library(laGP)
library(DiceKriging)
library(DiceOptim)
library(jsonlite)
library(glue)

goldstein.price <- function(X)
{
    if(is.null(nrow(X))) X <- matrix(X, nrow=1)
    m <- 8.6928
    s <- 2.4269
    x1 <- 4 * X[,1] - 2
    x2 <- 4 * X[,2] - 2
    a <- 1 + (x1 + x2 + 1)^2 * (19 - 14 * x1 + 3 * x1^2 - 14 *
        x2 + 6 * x1 * x2 + 3 * x2^2)
    b <- 30 + (2 * x1 - 3 * x2)^2 * (18 - 32 * x1 + 12 * x1^2 +
        48 * x2 - 36 * x1 * x2 + 27 * x2^2)
    f <- log(a * b)
    f <- (f - m)/s
    return(f)
}

toy.c1 <- function(X) {
  if (is.null(dim(X))) X <- matrix(X, nrow=1)
  c1 <- 3/2 - X[,1] - 2*X[,2] - 0.5*sin(2*pi*(X[,1]^2 - 2*X[,2]))
  return(cbind(c1, -apply(X, 1, branin) + 25))
}

parr <- function(X){
  if (is.null(dim(X))) X <- matrix(X, nrow=1)
  x1 <- (2 * X[,1] - 1)
  x2 <- (2 * X[,2] - 1)
  g <- (4-2.1*x1^2+1/3*x1^4)*x1^2 + x1*x2 + (-4+4*x2^2)*x2^2+3*sin(6*(1-x1)) + 3*sin(6*(1-x2))
  return(-g+6)
}

gsbp.constraints <- function(x){
  return(cbind(toy.c1(x), parr(x)-2))
}

## problem definition for AL
gsbpprob <- function(X, known.only=FALSE)
{
  if(is.null(nrow(X))) X <- matrix(X, nrow=1)
  if(known.only) stop("known.only not supported for this example")
  f <- goldstein.price(X)
  C <- gsbp.constraints(X)
  return(list(obj=f, c=cbind(C[,1], C[,2]/100, C[,3]/10))) # Dividing by these constants isn't mentioned in Picheny et al. (https://arxiv.org/pdf/1605.09466.pdf)
}

## set bounding rectangle for aquisitions
dim <- 2

B <- matrix(c(rep(0,dim),rep(1,dim)),ncol=2)

ncandf <- function(t) { 1000 }


for(x in 1:100) {
  ## run ALBO
  ALslack <- optim.auglag(gsbpprob, B, equal=c(0,1,1), fhat=TRUE, urate=5, slack=TRUE, ncandf=ncandf, start=10, end=150)
  write_json(ALslack, glue("results/gsbp/slack_al/data/run_{x}_results.json"), digits=NA)
}
```


