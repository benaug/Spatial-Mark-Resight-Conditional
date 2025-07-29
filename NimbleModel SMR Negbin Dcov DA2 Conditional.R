NimModel <- nimbleCode({
  #baseline density
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  D.beta1 ~ dnorm(0,sd=10)
  # D.beta0 ~ dnorm(0,sd=10)

  #detection function priors
  lam0 ~ dunif(0,15)
  # phi ~ dunif(0,100) #not really uninformative, puts a lot of prior mass on very low overdispersion
  phi ~ dgamma(1,0.1) #moderately informative for high to moderate overdispersion
  sigma ~ dunif(0,10)

  #sample type observation model priors (Dirichlet)
  alpha.marked[1] <- 1
  alpha.marked[2] <- 1
  alpha.marked[3] <- 1
  alpha.unmarked[1] <- 1
  alpha.unmarked[2] <- 1
  theta.marked[1:3] ~ ddirch(alpha.marked[1:3])
  theta.unmarked[1] <- 0
  theta.unmarked[2:3] ~ ddirch(alpha.unmarked[1:2])

  #Density model
  D.intercept <- D0*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells])
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  lambda.N <- D.intercept*pi.denom #Expected N
  N ~ dpois(lambda.N) #realized N in state space.
  
  #Marked individuals first
  for(i in 1:n.marked){
    #dunif() here implies uniform distribution within a grid cell
    #also tells nimble s's are in continuous space, not discrete
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    #categorical likelihood for this cell, equivalent to zero's trick
    #also disallowing s's in non-habitat
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    y.true[i,1:J] ~ dNBVector(mu=lam[i,1:J],phi=phi,K1D=K1D[1:J],z=z[i])
    #custom distribution that skips likelihood eval for the individuals currently with 0 captures.
    y.event[i,1:J,1:3] ~ dmulti2(y.true[i,1:J],prob=theta.marked[1:3],capcounts=capcounts[i])
  }#custom Metropolis-Hastings update for N.M/z[1:n.marked] 
  
  #Then unmarked individuals
  for(i in (n.marked+1):M){
    #dunif() here implies uniform distribution within a grid cell
    #also tells nimble s's are in continuous space, not discrete
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    #categorical likelihood for this cell, equivalent to zero's trick
    #also disallowing s's in non-habitat
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    y.true[i,1:J] ~ dNBVector(mu=lam[i,1:J],phi=phi,K1D=K1D[1:J],z=z[i])
    #custom distribution that skips likelihood eval for the individuals currently with 0 captures.
    y.event[i,1:J,2:3] ~ dmulti2(y.true[i,1:J],prob=theta.unmarked[2:3],capcounts=capcounts[i])
  }#custom Metropolis-Hastings update for N.UM/z[(n.marked+1):M]
  capcounts[1:M] <- Getcapcounts(ID=ID[1:n.samples],capcounts.ID=capcounts.ID[1:M])
  n.cap <- Getncap(capcounts=capcounts[1:M])

  #If you have telemetry
  for(i in 1:n.tel.inds){
    for(m in 1:n.locs.ind[i]){
      locs[i,m,1] ~ dnorm(s[tel.inds[i],1],sd=sigma)
      locs[i,m,2] ~ dnorm(s[tel.inds[i],2],sd=sigma)
    }
  }
})# end model