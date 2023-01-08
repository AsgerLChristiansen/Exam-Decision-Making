install.packages("pacman")
pacman::p_load(R2jags, parallel)

set.seed(1337)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#----------getting the data

#load control data
ctr_data <- read.table("216377/Module4/data/IGTdata_healthy_control.txt",header=TRUE)

#----------prepare data for jags models - want trial x subject arrays for choice and gain & loss ----
# identify and count unique subject IDs
subIDs <- unique(ctr_data$subjID)
nsubs <- length(subIDs)
ntrials_max <- 100

# all choices (x) and outcomes (X)
x_raw <- ctr_data$deck
X_raw <- ctr_data$gain + ctr_data$loss #note the sign!

#--- assign choices and outcomes in trial x sub matrix

#different number of trials across subjects. We'll need to fix this by padding arrays of < 100
#this is just so we can make the array
#then we'll also need to record number of valid trials for each sub, 
#then run the JAGS model on only valid trials

# empty arrays to fill
ntrials_all <- array(0,c(nsubs))
x_all <- array(0,c(nsubs,ntrials_max))
X_all <- array(0,c(nsubs,ntrials_max))

for (s in 1:nsubs) {
  
  #record n trials for subject s
  ntrials_all[s] <- length(x_raw[ctr_data$subjID==subIDs[s]])
  
  #pad trials with NA if n trials < maximum (i.e. 100)
  x_sub <- x_raw[ctr_data$subjID==subIDs[s]] 
  length(x_sub) <- ntrials_max
  
  X_sub <- X_raw[ctr_data$subjID==subIDs[s]] 
  length(X_sub) <- ntrials_max
  
  # assign arrays
  x_all[s,] <- x_sub
  X_all[s,] <- X_sub
  
}


#----------testing our data curation by running JAGS on one subject

# Now we'll fit one subject just to make sure everything works

x <- x_all[1,]
X <- X_all[1,]

ntrials <- ntrials_all[1]

# set up jags and run jags model on one subject
data <- list("x","X","ntrials") 
params<-c("w","A","theta","a","p")
temp_samples <- jags.parallel(data, inits=NULL, params,
                model.file ="216377/Module3/PVL.txt",
                n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=4)

# let's look at the posteriors for the parameters
par(mfrow=c(2,2))
plot(density(temp_samples$BUGSoutput$sims.list$w))
plot(density(temp_samples$BUGSoutput$sims.list$A))
plot(density(temp_samples$BUGSoutput$sims.list$theta))
plot(density(temp_samples$BUGSoutput$sims.list$a))

# Question: how would you expect the data to look on the basis of these posteriors?



###########################################################
#---------- run the hierarchical model on controls --------
###########################################################

x <- x_all
X <- X_all

ntrials <- ntrials_all

# set up jags and run jags model
data <- list("x","X","ntrials","nsubs") 
params<-c("mu_w","mu_A","mu_theta","mu_a","lambda_w","lambda_A","lambda_theta","lambda_a")

samples <- jags.parallel(data, inits=NULL, params,
                model.file ="216377/Module4/hier_PVL.txt",
                n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=4)

# let's look at the posteriors for the parameters
par(mfrow=c(2,2))
plot(density(samples$BUGSoutput$sims.list$mu_w))
plot(density(samples$BUGSoutput$sims.list$mu_A))
plot(density(samples$BUGSoutput$sims.list$mu_theta))
plot(density(samples$BUGSoutput$sims.list$mu_a))

# let's look at the posteriors for the parameters
par(mfrow=c(2,2))
plot(density(samples$BUGSoutput$sims.list$lambda_w))
plot(density(samples$BUGSoutput$sims.list$lambda_A))
plot(density(samples$BUGSoutput$sims.list$lambda_theta))
plot(density(samples$BUGSoutput$sims.list$lambda_a))

