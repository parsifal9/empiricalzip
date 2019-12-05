#' pEMAlgorithm
#'
#' Fits a Poisson by EM 
#' @param x The data.
#' @keywords cats
#' @export
#' @examples
#' pEMAlgorithm(x)

pEMAlgorithm <- function (x)
{
  K <- max(x)
  M <- length(x)
  
  pEstimate <- matrix(0, nrow=3, ncol=K)
  
  for(i in 1:K){
    C <- i
    count <- c(length(x[which(x==0)]), tabulate(x, nbins = C))
    n <- sum(count)
    lambdaEstimate <- max(0, mean(x[which(x!=0)]))
    epsilon <- 10^(-4)
    cf <- 10^(-10)
    hold.eta <- 0
    hold.lambda <- max(cf, lambdaEstimate)
    hold.theta <- 0
    A <- sum(c(0:C)*count) + f3(x=x, eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
    B <- n + f(eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
    lambdaEstimate <- A/(B + cf)
    
    while((abs(hold.lambda - lambdaEstimate)) > epsilon){
      hold.eta <- 0
      hold.lambda <- max(cf, lambdaEstimate)
      hold.theta <- 0
      A <- sum(c(0:C)*count) + f3(x=x, eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
      B <- n + f(eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
      lambdaEstimate <- A/(B + cf)
    }
    pEstimate[1,i] <- 0
    pEstimate[2,i] <- lambdaEstimate
    pEstimate[3,i] <- 0
  }
  class(pEstimate) <- "pEstimate"
  return(pEstimate)
}
