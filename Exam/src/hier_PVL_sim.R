hier_PVL_sim <- function(payoff,nsubs,ntrials,mu_w,mu_A,mu_a,mu_theta,
                         sigma_w, sigma_A, sigma_a, sigma_theta) {
  
  # arrays to populate for simulation
  x <- array(NA,c(nsubs,ntrials[1]))
  X <- array(NA,c(nsubs,ntrials[1]))
  Ev <- array(NA,c(nsubs,ntrials[1],4))
  
  for (s in 1:nsubs) {
    # free parameters - sampled from a normal distribution with group mean and sd
    # w <- rnorm(1,mu_w,sigma_w) # not yet truncating this
    # A <- rnorm(1,mu_A,sigma_A)
    # theta <- rnorm(1,mu_theta,sigma_theta)
    # a <- rnorm(1,mu_a,sigma_a)
    w <- rtruncnorm(1,0,,mu_w,sigma_w) # not yet truncating this
    A <- rtruncnorm(1,0,,mu_A,sigma_A)
    theta <- rtruncnorm(1,0,,mu_theta,sigma_theta)
    a <- rtruncnorm(1,0,,mu_a,sigma_a)
    
    # arrays to populate for simulation
    u <- array(NA,c(ntrials[s],4))
    Ev_update <- array(NA,c(ntrials[s],4))
    exp_p <- array(NA,c(ntrials[s],4))
    p <- array(NA,c(nsubs,ntrials[s],4))
    
    #--- plot prospect theory function
    #x <- seq(1,100,1)
    #y <- x^A
    #plot(x,y)
    
    x[s,1] <- rcat(1,c(.25,.25,.25,.25)) # assigning a "flat" probability structure to the first choice (i.e. random choice between the four decks)
    
    X[s,1] <- payoff[1, x[s, 1]] # assigning the payoff to first random choice
    
    Ev[s,1,] <- rep(0,4) # assigning zero as the expected value for all four decks at the first "random" choice 
    
    for (t in 2:ntrials[s]) {
      
      for (d in 1:4) {
        
        u[t,d] <- ifelse(X[s,t-1]<0,-w*abs(X[s,t-1])^A,X[s,t-1]^A)
        
        Ev_update[t,d] <- Ev[s,t-1,d] + (a * (u[t] - Ev[s,t-1,d]))
        
        Ev[s,t,d] <- ifelse(x[s,t-1]==d,Ev_update[t,d],Ev[s,t-1,d])
        
        exp_p[t,d] <- exp(theta*Ev[s,t,d])
        
      }
      
      for (d in 1:4) {
        p[s,t,d] <- exp_p[t,d]/sum(exp_p[t,])
      }
      
      x[s,t] <- rcat(1,p[s,t,])
      
      X[s,t] <- payoff[t,x[s, t]]
      
    }
  }
  
  result <- list(x=x,
                 X=X,
                 Ev=Ev)
  
  return(result)
  
  
  #turn back on when building
  #par(mfrow=c(2,2))
  #plot(Ev[,1])
  #plot(Ev[,2])
  #plot(Ev[,3])
  #plot(Ev[,4])
  #plot(x)
  
}
