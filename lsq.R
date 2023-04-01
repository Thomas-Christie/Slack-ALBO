library(laGP)
library(jsonlite)
library(glue)

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
  write_json(out, glue("results/lsq/slack_al/data/run_{x}_results.json"), digits=NA)
}


