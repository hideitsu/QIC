## Before performing experiment, you have to put data in ./Data/
## For reference, we only provide wine.data in this source code distribution

getData <- function(DataName="wine",seed=123,red=NULL){
  set.seed(seed)
  if(DataName=="wine"){
    tmp <- read.csv(file="./Data/wine.data",head=FALSE)
    classes <- as.numeric(tmp[,1])
    id <- list()
    for(i in 1:length(unique(classes))){
      tmpid <- which(classes==i);len <- length(tmpid)
      id[[i]] <- sample(tmpid, round(len*.9))
    }
    X <- as.matrix(tmp[unlist(id),-1]);colnames(X) <- rownames(X) <- NULL
    classes <- classes[unlist(id)]
  }
  
  return(list(x=X,classes=as.factor(classes)))
}


