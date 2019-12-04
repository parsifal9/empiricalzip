#' zipEMAlgorithm
#'
#' Fits a zero-inflated by EM 
#' @param x The data.
#' @keywords cats
#' @export
#' @examples
#' zipEMAlgorithm(x)

zipEMAlgorithm <- function (x)
{
  K <- max(x)
  M <- length(x)
  
  zipEstimate <- matrix(0, nrow=3, ncol=K)
  
  for(i in 1:K){
    C <- i
    count <- c(length(x[which(x==0)]), tabulate(x, nbins = C))
    n <- sum(count)
    etaEstimate <- max(0,mean((x==0)))
    lambdaEstimate <- max(0,mean(x)/(1 - mean((x==0))))
    epsilon <- 10^(-4)
    cf <- 10^(-10)
    hold.eta <- max(0,etaEstimate)
    hold.lambda <- max(cf,lambdaEstimate)
    hold.theta <- 0
    u1 <- count[1]*tau0(x=0, eta=hold.eta, lambda=hold.lambda, theta=hold.theta)
    u2 <- sum(count*tau1(x=c(0:C), eta=hold.eta, lambda=hold.lambda, theta=hold.theta))
    u3 <- f(eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
    etaEstimate <- u1/(u1 + u2 + u3 + cf)
    v1 <- sum(c(0:C)*count*tau1(x=c(0:C), eta=hold.eta, lambda=hold.lambda, theta=hold.theta))
    v2 <- f3(x=x, eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
    lambdaEstimate <- (v1 + v2)/(u2 + u3 + cf)
    
    while(((abs(hold.eta - etaEstimate)) > epsilon)||((abs(hold.lambda - lambdaEstimate)) > epsilon)){
      hold.eta <- max(0,etaEstimate)
      hold.lambda <- max(cf,lambdaEstimate)
      hold.theta <- 0
      u1 <- count[1]*tau0(x=0, eta=hold.eta, lambda=hold.lambda, theta=hold.theta)
      u2 <- sum(count*tau1(x=c(0:C), eta=hold.eta, lambda=hold.lambda, theta=hold.theta))
      u3 <- f(eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
      etaEstimate <- u1/(u1 + u2 + u3 + cf)
      v1 <- sum(c(0:C)*count*tau1(x=c(0:C), eta=hold.eta, lambda=hold.lambda, theta=hold.theta))
      v2 <- f3(x=x, eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
      lambdaEstimate <- (v1 + v2)/(u2 + u3 + cf)
    }
    zipEstimate[1,i] <- etaEstimate
    zipEstimate[2,i] <- lambdaEstimate
    zipEstimate[3,i] <- 0
  }
  return(zipEstimate)
}
