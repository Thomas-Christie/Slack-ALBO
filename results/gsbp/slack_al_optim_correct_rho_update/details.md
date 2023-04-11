Ran Slack-AL + Optim on GSBP problem with correct handling of penalty updates.

``` 
library(laGP)
library(DiceKriging)
library(DiceOptim)
library(jsonlite)
library(glue)
library(R.utils)
library(tgp)

new_rho_update <- function (obj, C, equal, init = 10, ethresh = 0.01)
{
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

new_auglag <- function(fn, B, fhat=FALSE, equal=FALSE, ethresh=1e-2, slack=FALSE,
  cknown=NULL, start=10, end=100, Xstart=NULL, sep=TRUE, ab=c(3/2,8),
  lambda=1, rho=NULL, urate=10, ncandf=function(t) { t }, dg.start=c(0.1,1e-6),
  dlim=sqrt(ncol(B))*c(1/100,10), Bscale=1, ey.tol=1e-2, N=1000,
  plotprog=FALSE, verb=2, ...)
{
  ## check start
  if(start >= end) stop("must have start < end")

  ## check sep and determine whether to use GP or GPsep commands
  if(sep) { newM <- newGPsep; mleM <- mleGPsep; updateM <- updateGPsep;
    alM <- laGP:::alGPsep; deleteM <- deleteGPsep; nd <- nrow(B) }
  else { newM <- newGP; mleM <- mleGP; updateM <- updateGP;
    alM <- laGP:::alGP; deleteM <- deleteGP; nd <- 1 }
  formals(newM)$dK <- TRUE;

  ## get initial designwhich.min
  X <- dopt.gp(start, Xcand=lhs(10*start, B))$XX
  ## X <- lhs(start, B)
  X <- rbind(Xstart, X)
  start <- nrow(X)

  ## first run to determine dimensionality of the constraint
  out <- fn(X[1,]*Bscale, ...)
  nc <- length(out$c)

  ## check equal argument
  if(length(equal) == 1) {
    if(!(equal %in% c(0,1)))
      stop("equal should be a logical scalar or vector of lenghth(fn(x)$c)")
    if(equal) equal <- rep(1, nc)
    else equal <- rep(0, nc)
  } else {
    if(length(equal) != nc)
      stop("equal should be a logical scalar or vector of lenghth(fn(x)$c)")
    if(! any(equal != 0 | equal != 1))
      stop("equal should be a logical scalar or vector of lenghth(fn(x)$c)")
  }
  equal <- as.logical(equal)

  ## allocate progress objects, and initialize
  prog <- obj <- rep(NA, start)
  C <- matrix(NA, nrow=start, ncol=nc)
  obj[1] <- out$obj; C[1,] <- out$c
  if(all(out$c[!equal] <= 0) && all(abs(out$c[equal]) < ethresh)) prog[1] <- out$obj
  else prog[1] <- Inf

  ## remainder of starting run
  for(t in 2:start) { ## now that fn is vectorized we can probably remove for
    out <- fn(X[t,]*Bscale, ...)
    obj[t] <- out$obj; C[t,] <- out$c
    ## update best so far
    if(all(out$c[!equal] <= 0) && all(abs(out$c[equal]) < ethresh) &&
      out$obj < prog[t-1]) prog[t] <- out$obj
    else prog[t] <- prog[t-1]
  }

  ## handle initial lambda and rho values
  if(length(lambda) == 1) lambda <- rep(lambda, nc)
  if(is.null(rho)) rho <- new_rho_update(obj, C, equal)
  else if(length(rho) != 1 || rho <= 0) stop("rho should be a positive scalar")

  ## calculate AL for data seen so far
  if(!slack) {  ## Original AL
    Ce <- matrix(rep(as.numeric(equal), start), nrow=start, byrow=TRUE)
    Ce[Ce != 0] <- -Inf
    Cm <- pmax(C, Ce)
  } else {  ## AL with slack variables
    S <- pmax(- C - rho*matrix(rep(lambda, start), nrow=start, byrow=TRUE), 0)
    S[, equal] <- 0
    Cm <- C+S
  }
  al <- obj + Cm %*% lambda + Cm^2 %*% rep(1/(2*rho), nc)

  ## best auglag seen so far
  ybest <- min(al)
  since <- 0

  ## best valid so far
  m2 <- prog[start]

  ## initializing constraint surrogates
  Cgpi <- rep(NA, nc)
  d <- matrix(NA, nrow=nc, ncol=nd)
  Cnorm <- rep(NA, nc)
  for(j in 1:nc) {
    if(j %in% cknown) { Cnorm[j] <- 1; Cgpi[j] <- -1 }
    else {
      Cnorm[j] <- max(abs(C[,j]))
      Cgpi[j] <- newM(X, C[,j]/Cnorm[j], dg.start[1], dg.start[2])
      d[j,] <- mleM(Cgpi[j], param="d", tmin=dlim[1], tmax=dlim[2], ab=ab,
                    verb=verb-1)$d
    }
  }
  ds <- matrix(rowMeans(d, na.rm=TRUE), nrow=1)

  ## possibly initialize objective surrogate
  if(fhat) {
    fnorm <- max(abs(obj))
    fgpi <- newM(X, obj/fnorm, dg.start[1], dg.start[2])
    df <- mleM(fgpi, param="d", tmin=dlim[1], tmax=dlim[2], ab=ab, verb=verb-1)$d
    dfs <- matrix(df, nrow=1)
  } else { fgpi <- -1; fnorm <- 1 }

  ## init for loop
  mei <- Inf
  new.params <- FALSE

  ## keeping track
  meis <- mei
  lambdas <- as.numeric(lambda)
  rhos <- rho

  ## iterating over the black box evaluations
  for(t in (start+1):end) {

    ## lambda and rho update
    if(!slack) {  ## Original AL
      Ce <- matrix(rep(as.numeric(equal), nrow(C)), nrow=nrow(C), byrow=TRUE)
      Ce[Ce != 0] <- -Inf
      Cm <- pmax(C, Ce)
      al <- obj + C %*% lambda + Cm^2 %*% rep(1/(2*rho), nc)

      if(!is.finite(m2)) { ## update if no valid has been found
        lambda.new <- rep(0, length=ncol(C))
        rho.new <- rho/2
        # rho.new <- auto.rho(obj, C, equal)
        # if(rho.new >= rho) rho.new <- 9*rho/10
        # else if(rho.new < rho/2) rho.new <- rho/2

      } else { ## update if any valid has been found
        ck <- C[which.min(al),]
        lambda.new <- pmax(0, lambda + (1/rho) * ck)
        if(any(ck[equal] > 0) || any(abs(ck[!equal]) > ethresh)) rho.new <- rho/2
        else rho.new <- rho
      }

    } else {  ## Slack Variable AL
      S <- pmax(- C - rho*matrix(rep(lambda, nrow(C)), nrow=nrow(C), byrow=TRUE), 0)
      S[, equal] <- 0
      Cm <- C+S
      al <- obj + Cm %*% lambda + Cm^2 %*% rep(1/(2*rho), nc)
      cmk <- Cm[which.min(al),]

      if(!is.finite(m2)) { ## update if no valid has been found
        ## same as in non-slack case above
        lambda.new <- rep(0, length=ncol(C))
        rho.new <- rho/2
        # ## lambda.new <- lambda + (1/rho) * cmk
        # rho.new <- auto.rho(obj, C, equal)
        # ## rho.new <- auto.rho(obj, Cm)
        # if(rho.new >= rho) rho.new <- 9*rho/10
        # else if(rho.new < rho/2) rho.new <- rho/2
      } else { ## update if any valid has been found
        ## specific to slack variable case
        lambda.new <- lambda + (1/rho) * cmk
        ## if(mean(cmk^2) > quantile(rowMeans(Cm^2), p=0.05)) rho.new <- rho/2
        if(any(cmk[!equal] > 0) || any(abs(cmk[equal]) > ethresh)) rho.new <- rho/2
        else rho.new <- rho
      }
    }

    ## printing progress
    if(any(lambda.new != lambda) || rho.new != rho) {
      if(verb > 0) {
        cat("updating La:")
        if(rho.new != rho) cat(" rho=", rho.new, sep="")
        if(any(lambda.new != lambda))
          cat(" lambda=(", paste(signif(lambda.new,3), collapse=", "),
            ")", sep="")
        cat("\n")
      }
      new.params <- TRUE
    } else new.params <- FALSE

    ## confirm update of augmented lagrangian
    lambda <- lambda.new; rho <- rho.new
    if(!slack) ## original AL
      ybest <- min(obj + C %*% lambda + Cm^2 %*% rep(1/(2*rho), nc)) ## orig
    else ybest <- min(obj + Cm %*% lambda + Cm^2 %*% rep(1/(2*rho), nc)) ## slack

    ## keep track of lambda and rho
    lambdas <- rbind(lambdas, lambda)
    rhos <- cbind(rhos, rho)

    ## rebuild surrogates periodically under new normalized responses
    if(t > (start+1) && (t %% urate == 0)) {

      ## constraint surrogates
      Cnorm <- apply(abs(C), 2, max)
      print(glue("New Norm: {Cnorm}"))
      for(j in 1:nc) {
        if(j %in% cknown) next;
        deleteM(Cgpi[j])
        d[j,d[j,] < dlim[1]] <- 10*dlim[1]
        d[j,d[j,] > dlim[2]] <- dlim[2]/10
        Cgpi[j] <- newM(X, C[,j]/Cnorm[j], d[j,], dg.start[2])
        d[j,] <- mleM(Cgpi[j], param="d", tmin=dlim[1], tmax=dlim[2],
                      ab=ab, verb=verb-1)$d
      }
      ds <- rbind(ds, rowMeans(d, na.rm=TRUE))

      ## possible objective surrogate
      if(fhat) {
        deleteM(fgpi)
        fnorm <- max(abs(obj))
        df[df < dlim[1]] <- 10*dlim[1]
        df[df > dlim[2]] <- dlim[2]/10
        fgpi <- newM(X, obj/fnorm, df, dg.start[2])
        df <- mleM(fgpi, param="d", tmin=dlim[1], tmax=dlim[2], ab=ab,
                   verb=verb-1)$d
        dfs <- rbind(dfs, df)
      } else { df <- NULL }

      new.params <- FALSE
    }

    ## random candidate grid
    ncand <- ncandf(t)
    if(!is.finite(m2) || fhat)
      XX <- lhs(ncand, B)
    else XX <- rbetter(ncand, B, sum(X[which(obj == m2)[1],]))
    ## NOTE: might be a version of rbetter for fhat

    ## calculate composite surrogate, and evaluate EI and/or EY
    eyei <- alM(XX, fgpi, fnorm, Cgpi, Cnorm, lambda, 1/(2*rho), ybest,
                slack, equal, N, fn, Bscale)
    eis <- eyei$ei; by <- "ei"
    mei <- max(eis)
    nzei <- sum(eis > 0)
    if(nzei <= ey.tol*ncand) { eis <- -(eyei$ey); by <- "ey"; mei <- Inf }
    else if(nzei <= 0.1*ncand) {
      if(!is.finite(m2) || fhat)
        XX <- rbind(XX, lhs(10*ncand, B))
      else XX <- rbind(XX, rbetter(10*ncand, B, sum(X[which(obj == m2)[1],])))
      eyei <- alM(XX[-(1:ncand),], fgpi, fnorm, Cgpi, Cnorm, lambda, 1/(2*rho), ybest,
                slack, equal, N, fn, Bscale)
      eis <- c(eis, eyei$ei)
      nzei <- nzei + sum(eis > 0)
      if(nzei <= ey.tol*ncand) stop("not sure")
      ncand <- ncand + 10*ncand
      mei <- max(eis)
    }
    meis <- c(meis, mei)

    ## plot progress
    if(!is.logical(plotprog) || plotprog) {
      par(mfrow=c(1,3+fhat))
      plot(prog, type="l", main="progress")
      if(is.logical(plotprog)) {
        if(length(eis) < 30) { span <- 0.5 } else { span <- 0.1 }
        g <- interp.loess(XX[,1], XX[,2], eis, span=span)
      } else g <- plotprog(XX[,1], XX[,2], eis)
      image(g, xlim=range(X[,1]), ylim=range(X[,2]), main="EI")
      valid <- apply(C, 1, function(x) { all(x <= 0) })
      if(is.matrix(valid)) valid <- apply(valid, 1, prod)
      points(X[1:start,1:2], col=valid[1:start]+3)
      points(X[-(1:start),1:2, drop=FALSE], col=valid[-(1:start)]+3, pch=19)
      matplot(ds, type="l", lty=1, main="constraint lengthscale")
      if(fhat) matplot(dfs, type="l", lty=1, main="objective lengthscale")
    }

    ## calculate next point
    m <- which.max(eis)
    xstar <- matrix(XX[m,], ncol=ncol(X))

    ## shift xstar via optim calculation
    if(slack > 1) {
      out <- laGP:::optim.alM(xstar, B, alM, by, fgpi, fnorm, Cgpi, Cnorm, lambda,
                       1/(2*rho), ybest, TRUE, equal, N, fn, Bscale)
      ## print(cbind(xstar, out$par))
      xstar <- out$par
    }

    ## progress meter
    if(verb > 0) {
      cat("t=", t, " ", sep="")
      cat(by, "=", eis[m]/Bscale, " (", nzei,  "/", ncand, ")", sep="")
      cat("; xbest=[", paste(signif(X[which(obj == m2)[1],],3), collapse=" "), sep="")
      cat("]; ybest (v=", m2, ", al=", ybest, ", since=", since, ")\n", sep="")
    }

    ## new run
    out <- fn(xstar*Bscale, ...)
    ystar <- out$obj; obj <- c(obj, ystar); C <- rbind(C, out$c)

    ## update GP fits
    X <- rbind(X, xstar)
    for(j in 1:nc) if(Cgpi[j] >= 0)
      updateM(Cgpi[j], xstar, out$c[j]/Cnorm[j], verb=verb-2)

    ## check if best valid has changed
    since <- since + 1
    if(all(out$c[!equal] <= 0) && all(abs(out$c[equal]) < ethresh) && ystar < prog[length(prog)]) {
      m2 <- ystar; since <- 0
    } ## otherwise m2 unchanged; should be the same as prog[length(prog)]
    prog <- c(prog, m2)

    ## check if best auglag has changed
    if(!slack) { ## original AL
      alstar <- out$obj + lambda %*% drop(out$c)
      ce <- as.numeric(equal); ce[ce != 0] <- -Inf
      alstar <- drop(alstar + rep(1/(2*rho),nc) %*% pmax(ce, drop(out$c))^2)
    } else {
      cs <-  pmax(-out$c - rho*lambda, 0)
      cs[equal] <- 0
      cps <- drop(out$c + cs)
      alstar <- out$obj + drop(lambda %*% cps + rep(1/(2*rho),nc) %*% cps^2)
    }
    if(alstar < ybest) { ybest <- alstar; since <- 0 }
  }

  ## delete GP surrogates
  for(j in 1:nc) if(Cgpi[j] > 0) deleteM(Cgpi[j])
  if(fhat) deleteM(fgpi)

  ## return output objects
  if(!fhat) df <- NULL
  return(list(prog=prog, mei=meis, obj=obj, X=X, C=C, d=d, df=df,
    lambda=as.matrix(lambdas), rho=as.numeric(rhos)))
}

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
  out <- new_auglag(gsbpprob, B, equal=c(0,1,1), fhat=TRUE, urate=5, slack=2, ncandf=ncandf, start=10, end=150)
  write_json(out, glue("results/gsbp/slack_al_optim_correct_rho_update/data/run_{x}_results.json"), digits=NA)
}
```