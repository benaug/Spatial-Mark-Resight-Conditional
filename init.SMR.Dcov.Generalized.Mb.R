e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.Dcov.Generalized.Mb <- function(data,inits=NA,M=NA,obsmod="poisson"){
  library(abind)
  #extract observed data
  y.mark <- data$y.mark
  y.mID <- data$y.mID #marked detections
  y.mnoID <- data$y.mnoID #marked with no ID samples
  y.um <- data$y.um #unmarked samples
  y.unk <- data$y.unk #unknown marked status samples
  n.marked <- data$n.marked
  X.mark <- as.matrix(data$X.mark)
  J.mark <- nrow(X.mark)
  K.mark <- data$K.mark
  K2D.mark <- data$K2D.mark
  X.sight <- as.matrix(data$X.sight)
  J.sight <- nrow(X.sight)
  K.sight <- data$K.sight
  K1D.sight <- data$K1D.sight
  locs <- data$locs
  
  xlim <- data$xlim
  ylim <- data$ylim
  
  ##pull out initial values
  p0.p <- inits$p0.p
  p0.c <- inits$p0.c
  lam0 <- inits$lam0
  sigma <- inits$sigma
  
  #augment marking process capture history
  y.aug <- array(0,dim=c(M,J.mark,K.mark))
  y.aug[1:n.marked,,] <- y.mark
  y.mark <- y.aug
  
  #restructure data into first and subsequent capture structures
  y.mark.p <- y.mark.c <- K1D.mark.p <- K1D.mark.c <- matrix(0,M,J.mark)
  for(i in 1:M){
    for(j in 1:J.mark){
      position <- Position(y.mark[i,j,],f=function(x){x>0})#occasion of first capture, NA if no capture
      if(is.na(position)){#no capture
        K1D.mark.p[i,j] <- sum(K2D.mark[j,1:K.mark]) #sum trap op
        K1D.mark.c[i,j] <- 0 #no subsequent capture events
      }else{#there is a first capture
        K1D.mark.p[i,j] <- sum(K2D.mark[j,1:position]) #sum trap op up to first capture event
        if(position<K.mark){#was first capture not last occasion?
          K1D.mark.c[i,j] <- sum(K2D.mark[j,(position+1):K.mark]) #sum trap op after first capture event for subsequent capture
        }else{#otherwise, no subsequent capture events
          K1D.mark.c[i,j] <- 0
        }
        y.mark.p[i,j] <- y.mark[i,j,position]
      }
    }
  }
  
  y.mark2D <- apply(y.mark,c(1,2),sum)
  y.mark.c <- y.mark2D-y.mark.p
  
  #I have not tested that the algorithm above is correct in all cases of trap operation entered.
  ##Check for bugs here##
  #I've never seen errors in this algorithm when the data is properly formatted. Errors may just be formatting errors.
  for(i in 1:M){
    for(j in 1:J.mark){
      y.check <- y.mark[i,j,]
      if(sum(y.check)>0){
        if(y.mark.p[i,j]!=1)stop("bug in y1")
        if(y.mark.c[i,j]!=(sum(y.check)-1))stop("bug in y2")
        first.cap.k.on <- which(y.check[K2D.mark[j,]==1]==1)[1] #first capture occasion counting only operable occasions
        if(K1D.mark.p[i,j]!=first.cap.k.on)stop("bug in K1D.p")
        if(K1D.mark.c[i,j]!=(sum(K2D.mark[j,1:K.mark])-first.cap.k.on))stop("bug in K1D.c")
      }else{
        if(y.mark.p[i,j]!=0)stop("bug in y1")
        if(y.mark.c[i,j]!=0)stop("bug in y2")
        if(K1D.mark.p[i,j]!=sum(K2D.mark[j,1:K.mark]))stop("bug in K1D.p")
        if(K1D.mark.c[i,j]!=0)stop("bug in K1D.c")
      }
    }
  }
  
  #assign random locations to assign latent ID samples to individuals
  s.init <- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2]))
  #but update s.inits for marked individuals before assigning latent detections
  y.both <- cbind(y.mark2D[1:n.marked,],y.mID)
  X.both <- rbind(X.mark,X.sight)
  idx <- which(rowSums(y.both)>0)
  for(i in idx){
    trps <- matrix(X.both[which(y.both[i,]>0),1:2],ncol=2,byrow=FALSE)
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
  #build plausible true sighting history to better initialize s
  D.sight <- e2dist(s.init, X.sight)
  lamd <- lam0*exp(-D.sight*D.sight/(2*sigma*sigma))
  y.true <- matrix(0,M,J.sight)
  y.true[1:n.marked,] <- y.mID
  y.event <- array(0,dim=c(M,J.sight,3))
  y.event[1:n.marked,,1] <- y.mID
  for(j in 1:J.sight){
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
  
  #update s.init given marking and sighting histories
  y.both <- cbind(y.mark2D,y.true)
  X.both <- rbind(X.mark,X.sight)
  idx <- which(rowSums(y.both)>0)
  for(i in idx){
    trps <- matrix(X.both[y.both[i,]>0,1:2],ncol=2,byrow=FALSE)
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
  
  #check starting logProbs
  D.mark <- e2dist(s.init, X.mark)
  pd.p <- p0.p*exp(-D.mark*D.mark/(2*sigma*sigma))
  pd.c <- p0.c*exp(-D.mark*D.mark/(2*sigma*sigma))
  D.sight <- e2dist(s.init, X.sight)
  lamd <- lam0*exp(-D.sight*D.sight/(2*sigma*sigma))
  
  #marking process
  logProb <- array(0,dim=c(M,J.mark))
  for(i in 1:M){
    for(j in 1:J.mark){
      logProb[i,j] <- dbinom(y.mark.p[i,j],K1D.mark.p[i,j],pd.p[i,j],log=TRUE)
      logProb[i,j] <- logProb[i,j] + dbinom(y.mark.c[i,j],K1D.mark.c[i,j],pd.c[i,j],log=TRUE)
    }
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite. Marking process.")
  
  #sighting process
  logProb <-  matrix(0,M,J.sight)
  if(obsmod=="poisson"){
    for(j in 1:J.sight){
      logProb[,j] <- dpois(y.true[,j],lamd[,j]*data$K1D.sight[j],log=TRUE)
    }
  }else if(obsmod=="negbin"){
    phi <- inits$phi
    for(j in 1:J.sight){
      p <- phi/(phi + lamd[,j])
      logProb[,j] <- dnbinom(y.true[,j],p=p,size=phi*data$K1D.sight[j],log = TRUE)
    }
  }else{
    print("obsmod must be 'poisson' or 'negbin'")
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood for sighting process not finite.")
  
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
    for(j in 1:J.sight){
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
    for(j in 1:J.sight){
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
    for(j in 1:J.sight){
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
  
  # y.event2 <- array(0,dim=c(M,J.sight,3))
  # y.event2[1:n.marked,,1] <- y.mID
  # for(l in 1:n.samples){
  #   y.event2[ID[l],this.j[l],event.type[l]] <- y.event2[ID[l],this.j[l],event.type[l]] + 1
  # }
  # all(y.event==y.event2)
  # y.true2 <- apply(y.event,c(1,2),sum)
  # all(y.true==y.true2)
  
  return(list(s=s.init,z=z.init,K1D.mark.p=K1D.mark.p,K1D.mark.c=K1D.mark.c,K1D.sight=K1D.sight,
              y.mark.p=y.mark.p,y.mark.c=y.mark.c,y.true=y.true,y.event=y.event,
              ID=ID,this.j=this.j,event.type=event.type,match=match,n.samples=n.samples,
              #y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,
              xlim=xlim,ylim=ylim,locs=locs,tel.inds=tel.inds,n.locs.ind=n.locs.ind))
}