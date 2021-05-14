## Evaluate clustering result with several indeces
source("DataHandle.r")
source("kmpp.r")
source("QIC.r")
library(clues)

DataName <- "wine"

alpha.cand <- seq(.05,.15,length=5)  ## canndidate values for alpha
anneal <- TRUE  ## if anneal <- TRUE, do not use cluster conditional entropy for evaluating the goodness of clustering to determine the best value of alpha

## maximum number of iteretion
maxit <- 10 ## sometimes better to use larger value

HAID <- NULL

for(s in 1:10){  ## 10-fold CV
  Dat <- getData(Data=DataName,seed=s)

  X <- Dat$x
  rmid <- which(apply(X,2,sd)==0)
  if(length(rmid)){
    X <- X[,-rmid]
  }
  X <- scale(X)
  
  K <- length(levels(Dat$classes))
  n <- length(Dat$classes)
  
  eps <- 1/(n*K)   ## convergence check
  p <- ncol(X); p <- min(p,30)  ## sometimes it is better to restrict dimensionality

  DM <- as.matrix(dist(X)); colnames(DM) <- NULL;rownames(DM) <- NULL
  ODM <- apply(DM,2,order)  ## ODM[-1,i] is the index of sorted data according to the distance from x_i
  SDM <- apply(DM,2,sort)  ## SDM[-1,i] is the sorted distance from x_i

  
  res <- QIC(DM,ODM,SDM,alpha.cand=alpha.cand,kmpp=TRUE,X=X,maxit=maxit,anneal=anneal,eps=eps)
  pred <- res$pred
  nW <- res$W

  HAID <- c(HAID,adjustedRand(as.numeric(Dat$classes),pred,randMethod="HA"))
}

## show average and standard deviation os HAID in 10 splits of training and test samples
print(c(mean(HAID),sd(HAID)))
