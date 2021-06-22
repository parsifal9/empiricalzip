#' empiricalfdr
#'
#' calculates the false discovery rate
#' @param estimate matrix with rows containing the estimates of eta, lambda, and theta
#' @param x   the count data in correct format
#' @param FDRLevel   the count data in correct format
#' @keywords cats
#' @export
#' @examples
#' empiricalfdr(estimate,x)

empiricalfdr<-function(estimate, x, FDRLevel=0.1){
    
    cf <- 10^(-10)
    count <- c(length(x[which(x==0)]), tabulate(x, nbins = max(x)))
    names(count) <- c(0:max(x))
    R <- rep(NA, 2)
    M<-length(x)  
    cut <- array(NA, max(x))
    pi <- array(NA, max(x))
    psi <- matrix(NA, nrow=2, ncol=max(x))
    xi <- array(NA, max(x))
    
    for(i in 1:max(x)){
        pi[i] <- min(1, sum(count[1:(i+1)])/(M*f0(estimate[1,i],estimate[2,i],estimate[3,i],C=i)))
        psi[1,i] <- length(x[which(x<=i)])*(-log(sum(count[1:(i+1)] + cf)/M)) + length(x[which(x>i)])*(-log(1 + cf - sum(count[1:(i+1)] + cf)/M)) # $\hat{C}_1$ Equation 6 
        psi[2,i] <- log(M)*length(x[which(x > i)]) - sum(count[(i + 2):(max(x) + 1)]*(log(count[(i + 2):(max(x) + 1)] + cf))) # $\hat{C}_2$ Equation 7 
        psi[2,max(x)] <- log(M)*length(x[which(x > i)]) - sum(count[max(x)]*(log(count[max(x)] + cf)))
        xi[i] <- length(x[which(x<=i)])*(-log(pi[i] + cf)) + length(x[which(x>i)])*(-log(1 + cf - pi[i]))
        #cut[i] <- psi[2,i] + negLogLh(x=x,estimate[1,i],estimate[2,i],estimate[3,i],C=i)   # $\hat{C}_1$ Equation 6
        cut[i] <- psi[2,i] + negLogLh(x=x,estimate[1,i],estimate[2,i],estimate[3,i],C=i)   # $\hat{C}_2$ Equation 7
        
    }
    
    C1 <- min(which(diff(sign(diff(cut)))>0)+2, which.min(cut))         # $\hat{C}_1$ Equation 6

    
    if(class(estimate)== "zigpEstimate"){
        AD <- min(max(x),ceiling(max((C1 + 1),estimate[2,C1]/(exp(estimate[3,C1] - 1) - estimate[3,C1]), log(length(x))/(estimate[3,C1] - 1 - log(estimate[3,C1])))))
        E<-rep(NA, 6)
        E[1] <- estimate[1,C1]
        E[2] <- estimate[2,C1]
        E[3] <- estimate[3,C1]
        E[4] <- pi[C1]
        E[5] <- C1
        E[6] <- AD
    } else if (class(estimate)== "zipEstimate"){
        AD <- min(max(x),ceiling(max((C1 + 1),estimate[2,C1], log(length(x)))))
        E<-rep(NA, 6)
        E[1] <- estimate[1,C1]
        E[2] <- estimate[2,C1]
        E[3] <- NA
        E[4] <- pi[C1]
        E[5] <- C1
        E[6] <- AD
    }  else if (class(estimate)== "gpEstimate"){
        AD <- min(max(x),ceiling(max((C1 + 1),estimate[2,C1]/(exp(estimate[3,C1] - 1) - estimate[3,C1]), log(length(x))/(estimate[3,C1] - 1 - log(estimate[3,C1])))))
        
        E<-rep(NA, 6)
        E[1] <- NA
        E[2] <- estimate[2,C1]
        E[3] <- estimate[3,C1]
        E[4] <- pi[C1]
        E[5] <- C1
        E[6] <- AD
    }else if (class(estimate)== "pEstimate"){
        AD <- min(max(x),ceiling(max((C1 + 1),estimate[2,C1], log(length(x)))))
        
        E<-rep(NA, 6)
        E[1] <- NA
        E[2] <- estimate[2,C1]
        E[3] <- NA
        E[4] <- pi[C1]
        E[5] <- C1
        E[6] <- AD
    }

    
    
    LFDR <- pi[C1]*gpdMixture(x=sort(x),estimate[1,C1],estimate[2,C1],estimate[3,C1])*f12(sort(x),length(x))
    sortx<-sort(x[which(x < (AD))])
    U1 <- pi[C1]*gpdMixture(x=sortx,estimate[1,C1],estimate[2,C1],estimate[3,C1])*f12(sortx,length(x))
    sortx<-sort(x[which(x  <=  C1)])
    U2 <- pi[C1]*gpdMixture(x=sortx,estimate[1,C1],estimate[2,C1],estimate[3,C1])*f12(sortx,length(x))
    R[1] <- sum(LFDR[!is.na(LFDR)] < FDRLevel)
    R[2] <- sum(U1[!is.na(U1)] < FDRLevel) - sum(U2[!is.na(U2)] < FDRLevel) + length(x[which(x > (AD-1))])
    BH <- RealMBHFDR(x,estimate[1,C1],estimate[2,C1],estimate[3,C1], FDRLevel = FDRLevel, M=length(x), piZero=pi[C1], frequency=count[which(count > 0)])
    Rejection <- cbind(R[1],R[2], BH)
    names(Rejection) <- c("One-stage procedure","Two-stage procedure","Storey s FDR")
    names(E) <- c("eta","lambda","theta","pi","C1","D(C1)")

    if (is.null(names(x))){
        names(x)<-c(1:length(x))
    }
    LFDR <- LFDR[match(names(x),names(x[order(x)]))]
    res<-list(LFDR, Rejection,E)
    names(res)<-c("LFDR", "Rejection","E")
    return(res)
}
