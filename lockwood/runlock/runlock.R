
## getting objective and constraints
## should separate into two functions: runlock and L
runlock <- function(rates=c(rep(1, 6)), known.only=FALSE)
{
  if(known.only) {
    if(is.null(nrow(rates))) rates <- matrix(rates, nrow=1)
    return(list(obj=rowSums(rates)))
  } else if(!is.null(nrow(rates)) && nrow(rates) != 1) stop("rates should be a vector") 
  write(c(6,rates * 10000),file="input.tst",ncol=1)
  system("./RunLock input.tst output.tst 1")
  output = scan(file="output.tst", quiet=TRUE)
  return(list(obj=output[1]/10000, c=output[-1]/1000))
}
