dCell <- nimbleFunction(
  run = function(x = double(0), pi.cell = double(0),log = integer(0)) {
    returnType(double(0))
    logProb <- log(pi.cell)
    return(logProb)
  }
)

#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(0)) {
    returnType(double(0))
    return(0)
  }
)

# Function to calculate detection rate, but skip when z=0
GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
     d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
     ans <- lam0*exp(-d2/(2*sigma^2))
     return(ans)
    }
  }
)
#Vectorized observation model that also prevents z from being turned off if an unmarked ind currently has samples.
#also skips likelihood eval when z=0
dNBVector <- nimbleFunction(
  run = function(x = double(1), mu = double(1), phi = double(0), K1D = double(1), z = double(0),
                 log = integer(0)){
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      p <- phi/(phi+mu)
      logProb <- sum(dnbinom(x,p=p,size=phi*K1D,log=TRUE))
      return(logProb)
    }
  }
)

#dummy random vector generator to make nimble happy
rNBVector <- nimbleFunction(
  run = function(n = integer(0),mu = double(1),phi = double(0), K1D = double(1), z = double(0)) {
    returnType(double(1))
    J <- nimDim(mu)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

#custom multinomial distribution to skip calcs when an ind has 0 samples (most of them!)
dmulti2 <- nimbleFunction(
  run = function(x = double(2), size = double(1), prob = double(1), capcounts = double(0),
                 log = integer(0)) {
    returnType(double(0))
    levels <- nimDim(prob)[1]
    J <- nimDim(size)[1]
    if(capcounts==0){
      return(0)
    }else{
      logProb <- 0
      for(j in 1:J){
        if(size[j]>0){
          logProb <- logProb + dmulti(x[j,1:levels], size=size[j], prob=prob, log = TRUE)
        }
      }
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rmulti2 <- nimbleFunction(
  run = function(n=integer(0),size = double(1), prob = double(1), capcounts = double(0)) {
    returnType(double(2))
    J <- nimDim(size)[1]
    out <- matrix(J,3,value=0)
    return(out)
  }
)

#calculates how many samples each individual is currently allocated.
Getcapcounts <- nimbleFunction(
  run = function(ID=double(1),capcounts.ID=double(1)){
    returnType(double(1))
    n.samples <- nimDim(ID)[1]
    capcounts <- capcounts.ID
    for(l in 1:n.samples){
      capcounts[ID[l]] <- capcounts[ID[l]] + 1
    }
    return(capcounts)
  }
)

#calculate number of captured individuals
Getncap <- nimbleFunction(
  run = function(capcounts=double(1)){
    returnType(double(0))
    M <- nimDim(capcounts)[1]
    nstate <- numeric(M, value = 0)
    for(i in 1:M){
      if(capcounts[i]>0){
        nstate[i] <- 1
      }
    }
    n.cap <- sum(nstate)
    return(n.cap)
  }
)

IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    n.marked <- control$n.marked
    M <- control$M
    J <- control$J
    event.type <- control$event.type
    match <- control$match
    n.samples <- control$n.samples
    this.j <- control$this.j
    calcNodes <- model$getDependencies(c("y.true","ID"))
  },
  run = function() {
    z <- model$z
    y.true <- model$y.true
    y.event <- model$y.event
    
    #precalculate log likelihoods.
    ll.y <- matrix(0,nrow=M,ncol=J)
    ll.y.event <- matrix(0,nrow=M,ncol=J)
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J) {
          #if phi is a function of individual or trap covariates, need to fix the lines below, e.g., size=model$phi[i,j]*model$K1D[j]
          p <- model$phi[1]/(model$phi[1] + model$lam[i,j])
          ll.y[i,j] <- dnbinom(y.true[i,j],p=p,size=model$phi[1]*model$K1D[j], log = TRUE)
        }
      }
    }
    for(i in 1:n.marked){
      if(z[i]==1){
        for(j in 1:J) {
          if(y.true[i,j]>0){
            ll.y.event[i,j] <- dmulti(y.event[i,j,1:3],y.true[i,j],model$theta.marked,log=TRUE)
          }
        }
      }
    }
    for(i in (n.marked+1):M){
      if(z[i]==1){
        for(j in 1:J) {
          if(y.true[i,j]>0){
            ll.y.event[i,j] <- dmulti(y.event[i,j,1:3],y.true[i,j],model$theta.unmarked,log=TRUE)
          }
        }
      }
    }
    ll.y.cand <- ll.y
    ll.y.event.cand <- ll.y.event
    ID.curr <- model$ID
    
    ###update IDs
    ID.cand <- ID.curr
    y.true.cand <- y.true
    y.event.cand <- y.event
    
    for(l in 1:n.samples){ #for all latent ID samples
      propprobs <- model$lam[1:M,this.j[l]]
      for(i in 1:M){ #zero out nonmatches (mnoID can't go to unmarked inds, um can't go to marked inds) and z=0
        if(!match[l,i] | z[i]==0){
          propprobs[i] <- 0
        }
      }
      propprobs <- propprobs/sum(propprobs)
      ID.cand[l] <- rcat(1,prob=propprobs)
      if(ID.cand[l]!=ID.curr[l]){ #abort if we propose the same individual, nothing changes
        swapped <- c(ID.curr[l],ID.cand[l])
        #new sample proposal probabilities
        forprob <- propprobs[swapped[2]]
        backprob <- propprobs[swapped[1]]
        #new y.true's - move sample from ID to ID.cand
        y.event.cand[ID.curr[l],this.j[l],event.type[l]] <- y.event[ID.curr[l],this.j[l],event.type[l]] - 1
        y.event.cand[ID.cand[l],this.j[l],event.type[l]] <- y.event[ID.cand[l],this.j[l],event.type[l]] + 1
        y.true.cand[ID.curr[l],this.j[l]] <- y.true[ID.curr[l],this.j[l]] - 1
        y.true.cand[ID.cand[l],this.j[l]] <- y.true[ID.cand[l],this.j[l]] + 1
        #if phi is a function of individual or trap covariates, need to fix these 2 lines below, e.g., size=model$phi[swapped[1],this.j[l]]*model$K1D[this.j[l]]
        p <- model$phi[1]/(model$phi[1]+model$lam[swapped[1],this.j[l]])
        ll.y.cand[swapped[1],this.j[l]] <- dnbinom(y.true.cand[swapped[1],this.j[l]],
                                                   p=p,size=model$phi[1]*model$K1D[this.j[l]],log=TRUE)
        p <- model$phi[1]/(model$phi[1]+model$lam[swapped[2],this.j[l]])
        ll.y.cand[swapped[2],this.j[l]] <- dnbinom(y.true.cand[swapped[2],this.j[l]],
                                                   p=p,size=model$phi[1]*model$K1D[this.j[l]],log=TRUE)
        #old ID event likelihood
        if(swapped[1]<=n.marked){#marked guy
          if(y.true.cand[swapped[1],this.j[l]]==0){
            ll.y.event.cand[swapped[1],this.j[l]] <- 0
          }else{
            ll.y.event.cand[swapped[1],this.j[l]] <- dmulti(y.event.cand[swapped[1],this.j[l],1:3],
                                                            y.true.cand[swapped[1],this.j[l]],model$theta.marked,log=TRUE)
          }
        }else{#unmarked guy
          if(y.true.cand[swapped[1],this.j[l]]==0){
            ll.y.event.cand[swapped[1],this.j[l]] <- 0
          }else{
            ll.y.event.cand[swapped[1],this.j[l]] <- dmulti(y.event.cand[swapped[1],this.j[l],1:3],
                                                            y.true.cand[swapped[1],this.j[l]],model$theta.unmarked,log=TRUE)
          }
        }
        #new ID event likelihood
        if(swapped[2]<=n.marked){#marked guy
          if(y.true.cand[swapped[2],this.j[l]]==0){
            ll.y.event.cand[swapped[2],this.j[l]] <- 0
          }else{
            ll.y.event.cand[swapped[2],this.j[l]] <- dmulti(y.event.cand[swapped[2],this.j[l],1:3],
                                                            y.true.cand[swapped[2],this.j[l]],model$theta.marked,log=TRUE)
          }
        }else{#unmarked guy
          if(y.true.cand[swapped[2],this.j[l]]==0){
            ll.y.event.cand[swapped[2],this.j[l]] <- 0
          }else{
            ll.y.event.cand[swapped[2],this.j[l]] <- dmulti(y.event.cand[swapped[2],this.j[l],1:3],
                                                            y.true.cand[swapped[2],this.j[l]],model$theta.unmarked,log=TRUE)
          }
        }
        #select sample to move proposal probabilities
        #P(select a sample of this type for this ID at this trap)
        focalprob <- y.event[swapped[1],this.j[l],event.type[l]]/n.samples
        focalbackprob <- y.event.cand[swapped[2],this.j[l],event.type[l]]/n.samples
        
        #sum log likelihoods and do MH step
        lp_initial <- sum(ll.y[swapped,this.j[l]]) + sum(ll.y.event[swapped,this.j[l]])
        lp_proposed <- sum(ll.y.cand[swapped,this.j[l]]) + sum(ll.y.event.cand[swapped,this.j[l]])
        log_MH_ratio <- (lp_proposed + log(backprob) + log(focalbackprob)) - (lp_initial + log(forprob) + log(focalprob))
        accept <- decide(log_MH_ratio)
        if(accept){
          y.event[swapped[1],this.j[l],event.type[l]] <- y.event.cand[swapped[1],this.j[l],event.type[l]]
          y.event[swapped[2],this.j[l],event.type[l]] <- y.event.cand[swapped[2],this.j[l],event.type[l]]
          y.true[swapped[1],this.j[l]] <- y.true.cand[swapped[1],this.j[l]]
          y.true[swapped[2],this.j[l]] <- y.true.cand[swapped[2],this.j[l]]
          ll.y[swapped[1],this.j[l]] <- ll.y.cand[swapped[1],this.j[l]]
          ll.y[swapped[2],this.j[l]] <- ll.y.cand[swapped[2],this.j[l]]
          ll.y.event[swapped[1],this.j[l]] <- ll.y.event.cand[swapped[1],this.j[l]]
          ll.y.event[swapped[2],this.j[l]] <- ll.y.event.cand[swapped[2],this.j[l]]
          ID.curr[l] <- ID.cand[l]
        }else{
          #set these back
          y.event.cand[swapped[1],this.j[l],event.type[l]] <- y.event[swapped[1],this.j[l],event.type[l]]
          y.event.cand[swapped[2],this.j[l],event.type[l]] <- y.event[swapped[2],this.j[l],event.type[l]]
          y.true.cand[swapped[1],this.j[l]] <- y.true[swapped[1],this.j[l]]
          y.true.cand[swapped[2],this.j[l]] <- y.true[swapped[2],this.j[l]]
          ID.cand[l] <- ID.curr[l]
        }
      }
    }
    #put everything back into the model$stuff
    model$y.true <<- y.true
    model$y.event <<- y.event
    model$ID <<- ID.curr
    model.lp.proposed <- model$calculate(calcNodes) #update logprob
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

#Required custom update for N/z
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    J <- control$J
    n.marked <- control$n.marked
    M <- control$M
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    lam.nodes <- control$lam.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function(){
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      if(updown==0){#subtract
        reject <- FALSE #we auto reject if you select a detected individual

        #find all z's currently on *excluding marked individuals*
        z.on <- which(model$z[(n.marked+1):M]==1) + n.marked
        n.z.on <- length(z.on)
        if(n.z.on>0){ #skip if no unmarked z's to turn off, otherwise nimble will crash
          pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
          pick <- z.on[pick]
          
          #prereject turning off marked individuals currently allocated samples
          if(model$capcounts[pick]>0){#is this an marked individual with samples?
            reject <- TRUE
          }
          if(!reject){
            #get initial logprobs for N and y
            lp.initial.N <- model$getLogProb(N.node)
            lp.initial.y <- model$getLogProb(y.nodes[pick])

            #propose new N/z
            model$N[1] <<-  model$N[1] - 1
            model$z[pick] <<- 0

            #turn off
            model$calculate(lam.nodes[pick])

            #get proposed logprobs for N and y
            lp.proposed.N <- model$calculate(N.node)
            lp.proposed.y <- model$calculate(y.nodes[pick])

            #MH step
            log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
            accept <- decide(log_MH_ratio)
            if(accept) {
              mvSaved["N",1][1] <<- model[["N"]]
              mvSaved["z",1][pick] <<- model[["z"]][pick]
              mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
            }else{
              model[["N"]] <<- mvSaved["N",1][1]
              model[["z"]][pick] <<- mvSaved["z",1][pick]
              model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
              model$calculate(y.nodes[pick])
              model$calculate(N.node)
            }
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M

          #find all z's currently off. Marked inds excluded here bc always on.
          z.off <- which(model$z[(n.marked+1):M]==0) + n.marked
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1

          #turn on
          model$calculate(lam.nodes[pick])

          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])

          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept){
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
            mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)