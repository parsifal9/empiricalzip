#' gpEMAlgorithm
#'
#' Fits a generalized Poisson by EM 
#' @param x The data.
#' @keywords cats
#' @export
#' @examples
#' gpEMAlgorithm(x)

gpEMAlgorithm <- function (x)
{
  K <- max(x)
  M <- length(x)
  
  gpEstimate <- matrix(0, nrow=3, ncol=K)
  
  for(i in 1:K){
    C <- i
    count <- c(length(x[which(x==0)]), tabulate(x, nbins = C))
    n <- sum(count)
    nz <- x[which(x!=0)]
    ex <- mean(nz[which(nz <= quantile(nz,0.5))])
    vx <- var(nz[which(nz <= quantile(nz,0.5))])
    phi <- sqrt(vx/(ex + 10^(-10)))
    thetaEstimate <- max(0, 1 - 1/phi)
    lambdaEstimate <- max(0, mean(x)/phi)
    epsilon <- 10^(-4)
    cf <- 10^(-10)
    hold.eta <- 0
    hold.lambda <- max(cf,lambdaEstimate)
    hold.theta <- max(0,thetaEstimate)
    u2 <- sum(count*tau1(x=c(0:C), eta=hold.eta, lambda=hold.lambda, theta=hold.theta))
    u3 <- f(eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
    u4 <- sum(c(0:C)*count*tau1(x=c(0:C), eta=hold.eta, lambda=hold.lambda, theta=hold.theta))
    u5 <- f3(x=x, eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
    A <- u4 + u5
    u6 <- sum(c(0, seq(0,C - 1,1))*count*tau2(x=c(0:C), lambda=hold.lambda, theta=hold.theta))
    u7 <- f5(x=x, eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
    u8 <- f7(x=x, eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
    B <- u2 + u3 + u6 + u7 - u8
    lambdaEstimate <- (B + cf)/(u2 + u3 + cf)
    v1 <- sum(c(0, seq(0,C - 1,1))*count*tau3(x=c(0:C), lambda=hold.lambda, theta=hold.theta))
    v2 <- f9(x=x, eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
    v3 <- f11(x=x, eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
    D <- v1 + v2 - v3
    thetaEstimate <- (D + cf)/(A + cf)
    
    while(((abs(hold.lambda-lambdaEstimate))>epsilon)||((abs(hold.theta-thetaEstimate))>epsilon)){
      hold.eta <- 0
      hold.lambda <- max(cf,lambdaEstimate)
      hold.theta <- max(0,thetaEstimate)
      u2 <- sum(count*tau1(x=c(0:C), eta=hold.eta, lambda=hold.lambda, theta=hold.theta))
      u3 <- f(eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
      u4 <- sum(c(0:C)*count*tau1(x=c(0:C), eta=hold.eta, lambda=hold.lambda, theta=hold.theta))
      u5 <- f3(x=x, eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
      A <- u4 + u5
      u6 <- sum(c(0, seq(0,C - 1,1))*count*tau2(x=c(0:C), lambda=hold.lambda, theta=hold.theta))
      u7 <- f5(x=x, eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
      u8 <- f7(x=x, eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
      B <- u2 + u3 + u6 + u7 - u8
      lambdaEstimate <- (B + cf)/(u2 + u3 + cf)
      v1 <- sum(c(0, seq(0,C - 1,1))*count*tau3(x=c(0:C), lambda=hold.lambda, theta=hold.theta))
      v2 <- f9(x=x, eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
      v3 <- f11(x=x, eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
      D <- v1 + v2 - v3
      thetaEstimate <- (D + cf)/(A + cf)
    }
    gpEstimate[1,i] <- 0
    gpEstimate[2,i] <- lambdaEstimate
    gpEstimate[3,i] <- thetaEstimate
  }
  return(gpEstimate)
}
