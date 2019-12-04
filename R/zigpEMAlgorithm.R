#' zigpEMAlgorithm
#'
#' Fits a zero-inflated generalized Poisson by EM 
#' @param x The data.
#' @keywords cats
#' @export
#' @examples
#' zigpEMAlgorithm(x)

zigpEMAlgorithm <- function (x)
{
  K <- max(x)
  M <- length(x)
  
  zigpEstimate <- matrix(0, nrow=3, ncol=K)
  
  for(i in 1:K){
    C <- i
    count <- c(length(x[which(x==0)]), tabulate(x, nbins = C))
    n <- sum(count)
    nz <- x[which(x!=0)]
    ex <- mean(x[which(x <= quantile(x,0.75))])
    vx <- var(x[which(x <= quantile(x,0.75))])
    er <- mean(x==0)/(1 - mean(x==0))
    phi <- sqrt(vx/ex - er*ex + 10^(-10))
    etaEstimate <- max(0,mean((x==0)))
    thetaEstimate <- min(max(0, 1 - 1/phi),1)
    lambdaEstimate <-  max(0, ex/((1 - mean(x==0))*phi))
    epsilon <- 10^(-4)
    cf <- 10^(-10)
    hold.eta <- max(0,etaEstimate)
    hold.lambda <- max(cf,lambdaEstimate)
    hold.theta <- max(cf,thetaEstimate)
    u1 <- count[1]*tau0(x=0, eta=hold.eta, lambda=hold.lambda, theta=hold.theta)
    u2 <- sum(count*tau1(x=c(0:C), eta=hold.eta, lambda=hold.lambda, theta=hold.theta))
    u3 <- f(eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
    etaEstimate <- u1/(u1 + u2 + u3 + cf)
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
    
    while(((abs(hold.eta-etaEstimate))>epsilon)||((abs(hold.lambda-lambdaEstimate))>epsilon)||((abs(hold.theta-thetaEstimate))>epsilon)){
      hold.eta <- max(0,etaEstimate)
      hold.lambda <- max(cf,lambdaEstimate)
      hold.theta <- max(cf,thetaEstimate)
      u1 <- count[1]*tau0(x=0, eta=hold.eta, lambda=hold.lambda, theta=hold.theta)
      u2 <- sum(count*tau1(x=c(0:C), eta=hold.eta, lambda=hold.lambda, theta=hold.theta))
      u3 <- f(eta=hold.eta, lambda=hold.lambda, theta=hold.theta, C=C, n=n)
      etaEstimate <- u1/(u1 + u2 + u3 + cf)
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
    zigpEstimate[1,i] <- etaEstimate
    zigpEstimate[2,i] <- lambdaEstimate
    zigpEstimate[3,i] <- thetaEstimate
  }
  return(zigpEstimate)
}
