# To run this script, you must specify either
# method <- "withscreening" # (this is algorithm 1)
# OR
# method <- "withoutscreening" # (this doesn't use the algorithms proposed in this paper)
# OR
# method <- "blockdiag" # (this is algorithm 2)
# OR
# method <- "MB" # (this is the approximation in Meinshausen and Buhlmann 2006)

# After running this script, be sure to completely exit R before you try to re-run it with a different choice of "method".
# Otherwise, your R will not properly load the other version of the R package since the 3 R packages' functions have the same names.

method <- "MB"
#method <- "withscreening"
#method <- "withoutscreening"
#method <- "blockdiag"
cat(method, fill=TRUE)

if(method=="MB" || method=="withoutscreening") library(glasso1.4)
if(method=="withscreening") library(glasso1.6)
if(method=="blockdiag") library(glasso1.7)
n <- 20
ps <- c(100, 800, 1500)
niter <- 10


FindLambdaWithGivenSparsity <- function(s, sparsity){
  ss <- s
  diag(ss) <- 0
  maxss <- apply(abs(ss), 2, max)
  lam <- quantile(maxss, sparsity)
  return(lam) # this is the lambda value such that the matrix will have sparsity% 0's.
}


msqrt <- function(Sigma){ # take a square root of a symmetric p.s.d. matrix
  if(sum((Sigma-t(Sigma))^2)>1e-7) stop("Sigma must be symmetric.")
  eig <- eigen(Sigma)
  return(eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors))
}

set.seed(10)


for(sparsitylevel in c(.2,.5,.9)){
  meantimes <- setimes <- NULL
  
  
  for(sim in 1:3){
    cat("Starting Simulation ", sim, "within sparsity level", sparsitylevel, fill=TRUE)
    if(method=="withscreening") filename <- paste("Sim",sim,"WITHSCREENING",sep="")
    if(method=="withoutscreening") filename <- paste("Sim",sim,"NOSCREENING",sep="")
    if(method=="blockdiag") filename <- paste("Sim",sim,"BLOCKDIAG",sep="")
    if(method=="MB") filename <- paste("Sim",sim,"MB",sep="")  
    times <- array(NA, dim=c(length(ps),  niter))
    
    for(iter in 1:niter){
      cat("Starting iteration ", iter, fill=TRUE)
      for(p in ps){
        cat("Starting p=", p, date(), fill=TRUE)
        x <- scale(matrix(rnorm(n*p), ncol=p), T, T)
        if(sim==2){
          Sigma <- matrix(0, nrow=p, ncol=p)
          Sigma[1:(p/2),1:(p/2)] <- 0.5
          diag(Sigma) <- 1
          Sigma.onehalf <- msqrt(Sigma)
          x <- x%*%Sigma.onehalf
        }
        if(sim==3){
          Sigma <- matrix(0.5, nrow=p, ncol=p)
          diag(Sigma) <- 1
          Sigma.onehalf <- msqrt(Sigma)
          x <- x%*%Sigma.onehalf
        }
        s <- t(x)%*%x/n
        lam <- FindLambdaWithGivenSparsity(s, sparsitylevel)
        if(method!="MB") times[ps==p, iter] <- system.time(out <- glasso(s, lam))[1]
        if(method=="MB") times[ps==p, iter] <- system.time(out <- glasso(s, lam, approx=TRUE))[1]
      }
    }
    meantimes <- cbind(meantimes, apply(times, 1,mean))
    setimes <- cbind(setimes, apply(times, 1,sd)/sqrt(niter))
  }
  if(method=="withscreening"){
    write.table(file=paste("SimTimesSCREENINGSparsityLevel",sparsitylevel,".txt",sep=""),
                matrix(paste(round(meantimes,3),"(",round(setimes,3),")",sep=""), ncol=3),row.names=ps, col.names=c("Sim1", "Sim2", "Sim3"), quote=FALSE)
  } else if(method=="withoutscreening"){
    write.table(file=paste("SimTimesNOSCREENINGSparsityLevel", sparsitylevel, ".txt",sep=""),
                matrix(paste(round(meantimes,3),"(",round(setimes,3),")",sep=""), ncol=3),row.names=ps, col.names=c("Sim1", "Sim2", "Sim3"), quote=FALSE)
  } else if(method=="blockdiag"){
    write.table(file=paste("SimTimesBLOCKDIAGSparsityLevel", sparsitylevel, ".txt",sep=""),
                matrix(paste(round(meantimes,3),"(",round(setimes,3),")",sep=""), ncol=3),row.names=ps, col.names=c("Sim1", "Sim2", "Sim3"), quote=FALSE)
  } else if(method=="MB"){
    write.table(file=paste("SimTimesMBSparsityLevel", sparsitylevel, ".txt",sep=""),
                matrix(paste(round(meantimes,3),"(",round(setimes,3),")",sep=""), ncol=3),row.names=ps, col.names=c("Sim1", "Sim2", "Sim3"), quote=FALSE)
  }
  
}

rm(list=ls())

