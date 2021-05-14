## kmeans++ for initializing weights.
## the original kmpp code for R is borrowed from
## https://stat.ethz.ch/pipermail/r-help/2012-January/300051.html
## written by Hans Werner

library(pracma)

kmpp_QIC <- function(X, k){
  n <- nrow(X)
  C <- numeric(k)
  C[1] <- sample(1:n, 1)
  for (i in 2:k) {
    dm <- distmat(X, X[C, ]); 
    pr <- apply(dm, 1, min); pr[C] <- 0;
    pr[is.nan(pr)] <- 0;    pr[is.na(pr)] <- 0
    C[i] <- sample(1:n, 1, prob = pr)
  }

  W <- NULL
  for(i in 1:k){
    W <- rbind(W,exp(-0.1*sd(DM[C[i],])*DM[C[i],]))
  }
  return(W)
}

