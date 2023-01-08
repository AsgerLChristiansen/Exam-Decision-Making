install.packages("pacman")
pacman::p_load(pacman,R2jags, parallel, ggpubr, extraDistr, truncnorm)

set.seed(1337)
setwd("/work/Exam")

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#------ create task environment -------------------
# NB! mod(ntrials, nstruct) (aka. ntrials %% nstruct) must be 0
ntrials <- 100 # total number of trials in our payoff structure
nstruct <- 10 # size of our subdivisions for pseudorandomization
freq <- 0.5 # probability of our frequent losses (we have losses half of the time)
infreq <- 0.1 # probability of our infrequent losses (we have losses 1/10th of the time)
bad_r <- 100 # "bad" winnings
bad_freq_l <- -250 # "bad" frequent loss
bad_infreq_l <- -1250 # "bad" infrequent loss
good_r <- 50 # "good" winnings
good_freq_l <- -50 # "good" frequent loss
good_infreq_l <- -250 # "good" infrequent loss

# Bad frequent
A_R <- rep(bad_r, nstruct) # we win on every trials
A_L <- c(rep(bad_freq_l, nstruct*freq),rep(0,nstruct*(1-freq))) # we have losses half of the time

# Bad infrequent
B_R <- rep(bad_r, nstruct)
B_L <- c(rep(bad_infreq_l, nstruct*infreq),rep(0,nstruct*(1-infreq))) # we have losses 1/10th of the time

# Good frequent
C_R <- rep(good_r, nstruct)
C_L <- c(rep(good_freq_l, nstruct*freq),rep(0,nstruct*(1-freq)))

# Good infrequent
D_R <- rep(good_r, nstruct)
D_L <- c(rep(good_infreq_l, nstruct*infreq),rep(0,nstruct*(1-infreq)))

# create the pseudorandomized full payoff structure
A <- array(NA,ntrials) # setting up and empty array to be filled
B <- array(NA,ntrials)
C <- array(NA,ntrials)
D <- array(NA,ntrials)
for (i in 1:(ntrials/nstruct)) {
  A[(1+(i-1)*nstruct):(i*nstruct)] <- (A_R + sample(A_L)) # randomly shuffling the loss-array for every ten trials (and adding those losses to the winnings)
  B[(1+(i-1)*nstruct):(i*nstruct)] <- (B_R + sample(B_L))
  C[(1+(i-1)*nstruct):(i*nstruct)] <- (C_R + sample(C_L))
  D[(1+(i-1)*nstruct):(i*nstruct)] <- (D_R + sample(D_L))
}


payoff <- cbind(A,B,C,D)/100 # combining all four decks as columns with each 100 trials - dividing our payoffs by 100 to make the numbers a bit easier to work with

# let's look at the payoff
colSums(payoff) # the two bad decks should sum to -25 (i.e. -2500), and the two good ones to 25 (i.e. 2500)

###--------------Run full parameter recovery -------------
niterations <- 100
nsubs <- 38
ntrials_all <- rep(100, nsubs) # all 38 subs have 100 trials each

# mu
true_mu_w <- array(NA,c(niterations))
true_mu_A <- array(NA,c(niterations))
true_mu_theta <- array(NA,c(niterations))
true_mu_a <- array(NA,c(niterations))

infer_mu_w <- array(NA,c(niterations))
infer_mu_A <- array(NA,c(niterations))
infer_mu_theta <- array(NA,c(niterations))
infer_mu_a <- array(NA,c(niterations))

# sigma (SD for R) / lambda (precision for JAGS)
true_lambda_w <- array(NA,c(niterations))
true_lambda_A <- array(NA,c(niterations))
true_lambda_theta <- array(NA,c(niterations))
true_lambda_a <- array(NA,c(niterations))

infer_lambda_w <- array(NA,c(niterations))
infer_lambda_A <- array(NA,c(niterations))
infer_lambda_theta <- array(NA,c(niterations))
infer_lambda_a <- array(NA,c(niterations))
?dgamma
start_time = Sys.time()
for (i in 1:niterations) {
  ntrials <- ntrials_all
  
  # let's see how robust the model is. Does it recover all sorts of values?
  mu_w <- runif(1,.5,2.5)
  mu_A <- runif(1,0,1)
  mu_theta <- runif(1,0,2)
  mu_a <- runif(1,0,1)
  
  sigma_w <- runif(1,0,0.2)
  sigma_A <- runif(1,0,0.1)
  sigma_theta <- runif(1,0,0.2)
  sigma_a <- runif(1,0,0.1)
  
  source('hier_PVL_sim.R')
  PVL_sims <- hier_PVL_sim(payoff,nsubs,ntrials,mu_w,mu_A,mu_a,mu_theta,
                  sigma_w, sigma_A, sigma_a, sigma_theta)
  
  x <- PVL_sims$x
  X <- PVL_sims$X
  # set up jags and run jags model
  data <- list("x","X","ntrials","nsubs") 
  params<-c("mu_w","mu_A","mu_theta","mu_a","lambda_w","lambda_A","lambda_theta","lambda_a")
  samples <- jags.parallel(data, inits=NULL, params,
                           model.file ="hier_PVL.txt", n.chains=3, 
                           n.iter=3000, n.burnin=1000, n.thin=1, n.cluster=4)
  
  # mu
  true_mu_w[i] <- mu_w
  true_mu_A[i] <- mu_A
  true_mu_theta[i] <- mu_theta
  true_mu_a[i] <- mu_a
  
  # find maximum a posteriori
  Y <- samples$BUGSoutput$sims.list
  infer_mu_w[i] <- MPD(Y$mu_w)
  infer_mu_A[i] <- MPD(Y$mu_A)
  infer_mu_theta[i] <- MPD(Y$mu_theta)
  infer_mu_a[i] <- MPD(Y$mu_a)
  
  # sigma
  true_lambda_w[i] <- sigma_w
  true_lambda_A[i] <- sigma_A
  true_lambda_theta[i] <- sigma_theta
  true_lambda_a[i] <- sigma_a
  
  # find maximum a posteriori
  infer_lambda_w[i] <- MPD(Y$lambda_w)
  infer_lambda_A[i] <- MPD(Y$lambda_A)
  infer_lambda_theta[i] <- MPD(Y$lambda_theta)
  infer_lambda_a[i] <- MPD(Y$lambda_a)
  
  print(i)
  
}

end_time = Sys.time()
end_time - start_time


save(true_lambda_a, file = "parameter_recovery_1/true_lambda_a.rda")
save(true_lambda_A, file = "parameter_recovery_1/true_lambda_A.rda")
save(true_lambda_theta, file = "parameter_recovery_1/true_lambda_theta.rda")
save(true_lambda_w, file = "parameter_recovery_1/true_lambda_w.rda")
save(true_mu_a, file = "parameter_recovery_1/true_mu_a.rda")
save(true_mu_A, file = "parameter_recovery_1/true_mu_A.rda")
save(true_mu_theta, file = "parameter_recovery_1/true_mu_theta.rda")
save(true_mu_w, file = "parameter_recovery_1/true_mu_w.rda")




save(infer_lambda_a, file = "parameter_recovery_1/infer_lambda_a.rda")
save(infer_lambda_A, file = "parameter_recovery_1/infer_lambda_A.rda")
save(infer_lambda_theta, file = "parameter_recovery_1/infer_lambda_theta.rda")
save(infer_lambda_w, file = "parameter_recovery_1/infer_lambda_w.rda")
save(infer_mu_a, file = "parameter_recovery_1/infer_mu_a.rda")
save(infer_mu_A, file = "parameter_recovery_1/infer_mu_A.rda")
save(infer_mu_theta, file = "parameter_recovery_1/infer_mu_theta.rda")
save(infer_mu_w, file = "parameter_recovery_1/infer_mu_w.rda")





setwd("/work/Exam")
load("parameter_recovery_1/true_lambda_a.rda")
load("parameter_recovery_1/true_lambda_A.rda")
load("parameter_recovery_1/true_lambda_theta.rda")
load("parameter_recovery_1/true_lambda_w.rda")
load("parameter_recovery_1/true_mu_a.rda")
load("parameter_recovery_1/true_mu_A.rda")
load("parameter_recovery_1/true_mu_theta.rda")
load("parameter_recovery_1/true_mu_w.rda")

load("parameter_recovery_1/infer_lambda_a.rda")
load("parameter_recovery_1/infer_lambda_A.rda")
load("parameter_recovery_1/infer_lambda_theta.rda")
load("parameter_recovery_1/infer_lambda_w.rda")
load("parameter_recovery_1/infer_mu_a.rda")
load("parameter_recovery_1/infer_mu_A.rda")
load("parameter_recovery_1/infer_mu_theta.rda")
load("parameter_recovery_1/infer_mu_w.rda")





# let's look at some scatter plots
# plotting code courtesy of Lasse
source('parameter_recovery_1/recov_plot.R')
pl1 <- recov_plot(true_mu_w, infer_mu_w, c("true mu_w", "infer mu_w"), 'smoothed linear fit')
pl2 <- recov_plot(true_mu_A, infer_mu_A, c("true mu_A", "infer mu_A"), 'smoothed linear fit')
pl3 <- recov_plot(true_mu_theta, infer_mu_theta, c("true mu_theta", "infer mu_theta"), 'smoothed linear fit')
pl4 <- recov_plot(true_mu_a, infer_mu_a, c("true mu_a", "infer mu_a"), 'smoothed linear fit')
ggarrange(pl1, pl2, pl3, pl4)


# mu_w is over-estimated. mu_A is looking pretty well.

# mu_theta is under-estimated pretty severely. mu_a is looking reasonable.

# In short, we can't quite trust estimates of mu_w and mu_theta. This is problematic, but there's little to do
# about it in the short time left.

# sigma (aka. true_lambda) re-coded as precision
pl1 <- recov_plot( true_lambda_w,infer_lambda_w, c("infer lambda_w","true lambda_w"), 'smoothed linear fit')
pl2 <- recov_plot(true_lambda_A,infer_lambda_A,  c("infer lambda_A","true lambda_A"), 'smoothed linear fit')
pl3 <- recov_plot(true_lambda_theta,infer_lambda_theta, c("infer lambda_theta", "true lambda_theta"), 'smoothed linear fit')
pl4 <- recov_plot(true_lambda_a,infer_lambda_a,  c("infer lambda_a", "true lambda_a"), 'smoothed linear fit')
ggarrange(pl1, pl2, pl3, pl4)

# lambda (aka. infer_lambda) re-coded as SD
pl1 <- recov_plot(1/sqrt(infer_lambda_w), 1/sqrt(true_lambda_w), c("infer lambda_w","true lambda_w"), 'smoothed linear fit')
pl2 <- recov_plot(1/sqrt(infer_lambda_A), 1/sqrt(true_lambda_A), c("infer lambda_A","true lambda_A"), 'smoothed linear fit')
pl3 <- recov_plot(1/sqrt(infer_lambda_theta), 1/sqrt(true_lambda_theta), c("infer lambda_theta", "true lambda_theta"), 'smoothed linear fit')
pl4 <- recov_plot(1/sqrt(infer_lambda_a), 1/sqrt(true_lambda_a), c("infer lambda_a", "true lambda_a"), 'smoothed linear fit')
ggarrange(pl1, pl2, pl3, pl4)

# Note - ridiculous over- and underestimation of lambda values.

