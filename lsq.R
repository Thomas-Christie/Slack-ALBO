library(laGP)

## this example assumes a known linear objective; further examples
## are in the optim.auglag demo
## a test function returning linear objective evaluations and
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
## optimization (primarily) by EI, change 25 to 100 for
## 99% chance of finding the global optimum with value 0.6
if(require("interp")) { ## for plotprog=interp
  out <- optim.auglag(aimprob, B, end=25, plotprog=interp)
} else {
  out <- optim.auglag(aimprob, B, end=25)
}
## using the slack variable implementation which is a little slower
## but more precise; slack=2 augments random search with L-BFGS-B
out2 <- optim.auglag(aimprob, B, end=25, slack=TRUE)
## Not run:
out3 <- optim.auglag(aimprob, B, end=25, slack=2)
## End(Not run)
