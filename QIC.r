## Quantile-based Information theoretic Clustering (QIC)
## 
## coded by Hideitsu Hino, 2014/01/02
## 
library(MASS)


weps <- 10e-10 ## lower bound of weight


## entropy estimator (up to constant)
SH <- function(DM, W=NULL, seed=123, sigma=NULL, alpha=.1,ID=NULL){
  W <- W/sum(W)
  id <- which(W < weps)
  diag(DM) <- 1

  if(length(id)>0){
    DM <- DM[-id,-id]
    W <- W[-id]
  }
  if(is.matrix(DM)){
    N <- nrow(DM)  ## data size
  }else{
    N <- length(DM)
  }

  ## weighted mean quantile
  if(N==1){
    print("Can not estimate entropy with only one datum.")
    return(NULL)
  }
  id <- which(DM==0,arr=T)[,1]
  if(length(id)!=0){
    DM <- DM[-id,-id]; W <- W[-id]
  }
  reg.w <- W/(1-W)
  H <- as.numeric(t(reg.w) %*% log(DM) %*% W)
  return(H)
}


### Main function
## DM: distance matrix
## ODM: ordered distance matrix
## SDM: sorted distance matrix
## maxit: maximum number of iteration
## alpha.cand: candidate values of alpha
## kmpp: if TRUE, k-means++ is used for initialization
## X: original data matrix (used for k-means++)
## anneal: if TRUE, alpha is gradually decreased. Else, optimal alpha is found by estimating cluster conditional entropy.
QIC <- function(DM,ODM,SDM,maxit=10,alpha.cand=0.1,kmpp=TRUE,X,anneal=TRUE,eps=1e-5){
  Wlist <- list();bestH <- Inf;bestid <- 1;predList <- list()

  ## constant term of the entropy estimator
  appC <- function(x){
    ((pi*exp(1)/(x/2))^(x/2))/sqrt(pi*x) - log(alpha)
  }
  
  if(anneal){
    alpha <- .9
    maxit <- 25
    if(kmpp){
      nW <- kmpp_QIC(X,K)
    }else{
      nW <- matrix(runif(K*n),ncol=n,nrow=K);nW <- nW %*% diag(1/colSums(nW))
    }
    ## find the alpha-quantile point. Can be improved...
    recalcLD <- function(W){
      tildeW <- diag(1/rowSums(W))%*%W
      alphaDM <- NULL
      for(k in 1:K){
        alphaD <- NULL
        for(i in 1:n){
          id <- which(((cumsum(tildeW[k,ODM[-1,i]]) >= alpha)==TRUE))[1]
          if(is.na(id)){
            id <- ODM[-1,i][length(ODM[-1,i])]
          }
          alphaD <- c(alphaD,SDM[-1,i][id])
        }
        alphaDM <- rbind(alphaDM, (alphaD))
      }
      return(alphaDM)
    }

    ## optimize the weight matrix
    nWold <- nW%*%diag(1/colSums(nW))
    prev.dif <- n*K
    for(i in 1:maxit){
      alpha <- max(0.01, alpha*0.8)
      alphaDM <- recalcLD(nW)
      if(sum(is.na(alphaDM))){break}
      tmp <- alphaDM^(-p)
      if(sum(is.infinite(tmp))){
        tmp <- (alphaDM+10^(-10))^(-1)
      }
      nW <- tmp%*%diag(1/colSums(tmp))
      pred <- apply(nW,2,which.max)
      current.dif <- mean(abs(nWold-nW))
      if( current.dif < eps){
        break
      }
      nWold <- nW
      prev.dif <- current.dif
    }
    
    bestid <- 1
    Wlist[[length(Wlist)+1]] <- nW
    predList[[length(predList)+1]] <- pred
    
    return(list(pred=predList[[bestid]],W=Wlist[[bestid]],best.alpha=alpha.cand[bestid]))
    
  }else{   ## when using cluster conditional entropy for evaluating goodness of cluster assignments
    for(r in 1:length(alpha.cand)){
      alpha <- alpha.cand[r]
      if(kmpp){
        nW <- kmpp_QIC(X,K)
      }else{
        nW <- matrix(runif(K*n),ncol=n,nrow=K);nW <- nW %*% diag(1/colSums(nW))
      }
      ## find the alpha-quantile point. Can be improved...
      recalcLD <- function(W){
        tildeW <- diag(1/rowSums(W))%*%W
        alphaDM <- NULL
        for(k in 1:K){
          alphaD <- NULL
          for(i in 1:n){
            id <- which(((cumsum(tildeW[k,ODM[-1,i]]) >= alpha)==TRUE))[1]
            if(is.na(id)){
              id <- ODM[-1,i][length(ODM[-1,i])]
            }
            alphaD <- c(alphaD,SDM[-1,i][id])
          }
          alphaDM <- rbind(alphaDM, (alphaD))
        }
        return(alphaDM)
      }
      


      ## optimize the weight matrix
      nWold <- nW%*%diag(1/colSums(nW))
      prev.dif <- n*K
      for(i in 1:maxit){
        alphaDM <- recalcLD(nW)
        if(sum(is.na(alphaDM))){break}
        tmp <- alphaDM^(-p)
        if(sum(is.infinite(tmp))){
          tmp <- (alphaDM+10^(-10))^(-1)
        }
        nW <- tmp%*%diag(1/colSums(tmp))
        pred <- apply(nW,2,which.max)
        current.dif <- mean(abs(nWold-nW))
        if( current.dif < eps){
          break
        }
        nWold <- nW
        prev.dif <- current.dif
      }
      
      if(length(alpha.cand)==1){
        bestid <- 1
        
        Wlist[[length(Wlist)+1]] <- nW
        predList[[length(predList)+1]] <- pred
      }else{
        H <- 0
        for(k in 1:K){
          tmp <- SH(DM,nW[k,]) * sum(nW[k,])/n
          H <- H+ifelse(length(tmp)==0,100,tmp)
        }
        
        Wlist[[length(Wlist)+1]] <- nW
        predList[[length(predList)+1]] <- pred
        
        if(bestH > H){
          bestH <- H
          bestid <- r
        }
      }
    }
    return(list(pred=predList[[bestid]],W=Wlist[[bestid]],best.alpha=alpha.cand[bestid]))
  }
}
