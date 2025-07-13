e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.Dcov <- function(data,inits=NA,M=NA,obsmod="poisson"){
  library(abind)
  #extract observed data
  y.mID <- data$y.mID #marked detections
  y.mnoID <- data$y.mnoID #marked with no ID samples
  y.um <- data$y.um #unmarked samples
  y.unk <- data$y.unk #unknown marked status samples
  n.marked <- data$n.marked
  X <- as.matrix(data$X)
  J <- nrow(X)
  K <- data$K
  K1D <- data$K1D
  locs <- data$locs
  
  xlim <- data$xlim
  ylim <- data$ylim
  
  ##pull out initial values
  lam0 <- inits$lam0
  sigma <- inits$sigma
  theta.marked <- inits$theta.marked
  theta.unmarked <- inits$theta.unmarked
 
  #assign random locations to assign latent ID samples to individuals
  s.init <- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2]))
  #but update s.inits for marked individuals before assigning latent detections
  idx <- which(rowSums(y.mID)>0)
  for(i in idx){
    trps <- matrix(X[which(y.mID[i,]>0),1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s.init[i,] <- trps
    }
  }
  #update using telemetry if you have it
  if(!is.null(dim(data$locs))){
    max.locs <- dim(locs)[2]
    if(n.marked>1){
      tel.inds <- which(rowSums(is.na(locs[,,1]))<max.locs)
      n.locs.ind <- rowSums(!is.na(locs[,,1]))
    }else{
      tel.inds <- which(sum(is.na(locs[,,1]))<max.locs)
      n.locs.ind <- sum(!is.na(locs[,,1]))
    }
    print("using telemetry to initialize telemetered s. Remove from data if not using in the model.")
    #update s starts for telemetry guys
    for(i in tel.inds){
      if(n.locs.ind[i]>1){
        s.init[i,] <- colMeans(locs[i,1:n.locs.ind[i],])
      }else{
        s.init[i,] <- locs[i,1,]
      }
      #make sure s is in state space
      if(s.init[i,1]<xlim[1]){
        s.init[i,1] <- xlim[1] + 0.01
      }
      if(s.init[i,1]>xlim[2]){
        s.init[i,1] <- xlim[2] - 0.01
      }
      if(s.init[i,2]<ylim[1]){
        s.init[i,2] <- ylim[1] + 0.01
      }
      if(s.init[i,2]>ylim[2]){
        s.init[i,2] <- ylim[2] - 0.01
      }
    }
    n.locs.ind <- n.locs.ind[tel.inds]
  }else{
    tel.inds <- NA
    n.locs.ind <- NA
  }
  
  D <- e2dist(s.init, X)
  lamd <- lam0*exp(-D*D/(2*sigma*sigma))
  y.true <- matrix(0,M,J)
  y.true[1:n.marked,] <- y.mID
  y.event <- array(0,dim=c(M,J,3))
  y.event[1:n.marked,,1] <- y.mID
  for(j in 1:J){
    #add marked no ID
    prob <- lamd[1:n.marked,j]
    prob <- prob/sum(prob)
    add <- rmultinom(1,y.mnoID[j],prob=prob)
    y.true[1:n.marked,j] <- y.true[1:n.marked,j] + add
    y.event[1:n.marked,j,2] <- add
    #add unmarked
    prob <- c(rep(0,n.marked),lamd[(n.marked+1):M,j])
    prob <- prob/sum(prob)
    add <- rmultinom(1,y.um[j],prob=prob)
    y.true[,j] <- y.true[,j] + add
    y.event[,j,2] <- y.event[,j,2] + add
    #add unk
    prob <- lamd[,j]
    prob <- prob/sum(prob)
    add <- rmultinom(1,y.unk[j],prob=prob)
    y.true[,j] <- y.true[,j] + add
    y.event[,j,3] <- y.event[,j,3] + add
  }
  
  z.init <- 1*(rowSums(y.true)>0)
  z.init[1:n.marked] <- 1
  
  #update s for individuals assigned samples, skip marked guys
  y.true2D <- y.true
  idx <- which(rowSums(y.true2D[(n.marked+1):M,])>0) + n.marked
  for(i in idx){
    trps <- matrix(X[y.true2D[i,]>0,1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s.init[i,] <- trps
    }
  }
  
  #If using a habitat mask, move any s's initialized in non-habitat above to closest habitat
  e2dist  <-  function (x, y){
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
  }
  getCell  <-  function(s,res,cells){
    cells[trunc(s[1]/res)+1,trunc(s[2]/res)+1]
  }
  alldists <- e2dist(s.init,data$dSS)
  alldists[,data$InSS==0] <- Inf
  for(i in 1:M){
    this.cell <- data$cells[trunc(s.init[i,1]/data$res)+1,trunc(s.init[i,2]/data$res)+1]
    if(data$InSS[this.cell]==0){
      cands <- alldists[i,]
      new.cell <- which(alldists[i,]==min(alldists[i,]))
      s.init[i,] <- data$dSS[new.cell,]
    }
  }
  
  D <- e2dist(s.init, X)
  lamd <- lam0*exp(-D*D/(2*sigma*sigma))
  
  
  #check starting logProbs
  logProb <-  matrix(0,M,J)
  if(obsmod=="poisson"){
    for(j in 1:J){
      logProb[,j] <- dpois(y.true[,j],lamd[,j]*data$K1D[j],log=TRUE)
    }
  }else if(obsmod=="negbin"){
    phi <- inits$phi
    for(j in 1:J){
      p <- phi/(phi + lamd[,j])
      logProb[,j] <- dnbinom(y.true[,j],p=p,size=phi*data$K1D[j],log = TRUE)
    }
  }else{
    print("obsmod must be 'poisson' or 'negbin'")
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite.")
  #reformat for conditional sampler
  n.samples <- sum(y.mnoID) + sum(y.um) + sum(y.unk)
  event.type <- c(rep("markednoID",sum(y.mnoID)),
                 rep("um",sum(y.um)),
                 rep("unk",sum(y.unk)))
  match <- matrix(FALSE,nrow=n.samples,ncol=M)
  #marked no ID samples can only match marked guys
  idx <- which(event.type=="markednoID")
  match[idx,1:n.marked] <- TRUE
  #unmarked samples can only match unmarked guys
  idx <- which(event.type=="um")
  match[idx,(n.marked+1):M] <- TRUE
  #unk samps can match anyone
  idx <- which(event.type=="unk")
  match[idx,] <- TRUE
  ID <- this.j <- rep(NA,n.samples)
  #mnoID
  idx <- 1
  for(i in 1:n.marked){
    for(j in 1:J){
      if(y.event[i,j,2]>0){
        for(l in 1:y.event[i,j,2]){
          ID[idx] <- i
          this.j[idx] <- j
          idx <- idx + 1
        }
      }
    }
  }
  #um
  for(i in (n.marked+1):M){
    for(j in 1:J){
      if(y.event[i,j,2]>0){
        for(l in 1:y.event[i,j,2]){
          ID[idx] <- i
          this.j[idx] <- j
          idx <- idx + 1
        }
      }
    }
  }
  #unk
  for(i in 1:M){
    for(j in 1:J){
      if(y.event[i,j,3]>0){
        for(l in 1:y.event[i,j,3]){
          ID[idx] <- i
          this.j[idx] <- j
          idx <- idx + 1
        }
      }
    }
  }
  
  #redefining event.type here to use in ID update
  event.type <- c(rep(2,sum(y.mnoID)),
                 rep(2,sum(y.um)),
                 rep(3,sum(y.unk)))
  
  # y.event2 <- array(0,dim=c(M,J,3))
  # y.event2[1:n.marked,,1] <- y.mID
  # for(l in 1:n.samples){
  #   y.event2[ID[l],this.j[l],event.type[l]] <- y.event2[ID[l],this.j[l],event.type[l]] + 1
  # }
  # all(y.event==y.event2)
  # y.true2 <- apply(y.event,c(1,2),sum)
  # all(y.true==y.true2)
  
  return(list(s=s.init,z=z.init,K1D=K1D,
              y.true=y.true,y.event=y.event,
              ID=ID,this.j=this.j,event.type=event.type,match=match,n.samples=n.samples,
              # y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,
              xlim=xlim,ylim=ylim,locs=locs,tel.inds=tel.inds,n.locs.ind=n.locs.ind))

}