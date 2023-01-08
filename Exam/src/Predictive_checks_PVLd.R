install.packages("pacman")
pacman::p_load(R2jags, parallel, ggplot2, tidyverse)

setwd("/work/Exam")
set.seed(1337)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#----------getting the data
#load control data from Maia & McClelland (undergraduate students used as controls)
ctr_choices <- read.table("raw_data/choice_100.txt",header=TRUE)
ctr_wins <- read.table("raw_data/wi_100.txt", header = TRUE)
ctr_losses <- read.table("raw_data/lo_100.txt", header = TRUE)
ctr_outcomes = ctr_wins + ctr_losses

ctr_authors <- read.table("raw_data/index_100.txt", header = TRUE)
ctr_maia <- ctr_authors %>% filter(Study == "Maia")
ctr_IDs = c(ctr_maia$Subj)

ctr_choices_maia = ctr_choices[ctr_IDs,]
ctr_outcomes_maia = ctr_outcomes[ctr_IDs,]

# Check if there are any NA's
NA_check = c(is.na(ctr_choices_maia[,100]))
print(NA_check) # Nope, all good

# Turn into data variables
x_all = ctr_choices_maia
X_all = ctr_outcomes_maia

nsubs <- as.numeric(length(x_all[,1]))
ntrials_all <- rep(100, nsubs) # all 38 subs have 100 trials each

ntrials <- ntrials_all[1]


# Scaling the payoffs (cuz we do so for the ORL to help better estimate the learning parameter)
X_all <- X_all/100


#----------testing our data curation by running JAGS on one subject

# Now we'll fit three subjects just to make sure everything works

x <- x_all[1,] # selecting choices from subj 1
X <- X_all[1,] # selecting payoffs from subj 1

nsubs <- 1
# set up jags and run jags model on one subject

# set up jags and run jags model
data <- list("x","X","ntrials","nsubs") 
params<-c("mu_w","mu_A","mu_theta","mu_a","lambda_w","lambda_A","lambda_theta","lambda_a", "p")
temp_samples <- jags.parallel(data, inits=NULL, params,
                         model.file ="parameter_recovery_2/hier_PVL_2.txt", n.chains=3, 
                         n.iter=3000, n.burnin=1000, n.thin=1, n.cluster=4)


# let's look at the posteriors for the parameters
par(mfrow=c(2,2))
plot(density(temp_samples$BUGSoutput$sims.list$mu_w))
plot(density(temp_samples$BUGSoutput$sims.list$mu_A))
plot(density(temp_samples$BUGSoutput$sims.list$mu_theta))
plot(density(temp_samples$BUGSoutput$sims.list$mu_a))

# Question: how would you expect the data to look on the basis of these posteriors?


#----------Posterior predictive checks of descriptive accuracy

# Posterior prediction - start by looking at posteriors for mu_w parameter

p_post <- temp_samples$BUGSoutput$sims.list$p # probabilities as the outcome from softmax


#plot probability of each deck on trial 84 for participant 2
par(mfrow=c(2,2))
plot(density(p_post[,1,84,1]))
plot(density(p_post[,1,84,2]))
plot(density(p_post[,1,84,3]))
plot(density(p_post[,1,84,4]))

# It would appear this participant has a very strong preference for deck 3 at trial 84.

# which option will be chosen?
x[1, 85] # 1.
# is this a good prediction?
# Well, no.


# let's write a loop that loop and see how the model goes at predicting responses for all trials 
x_predict <- array(ntrials-1) # We can't predict the first trial.

for (t in 1:(ntrials-1)){
  
  p_predict <- c(
    MPD(p_post[,1,t,1]),
    MPD(p_post[,1,t,2]),
    MPD(p_post[,1,t,3]),
    MPD(p_post[,1,t,4])
  )

  x_predict[t] <- which.max(p_predict)
}

x_compare = as.numeric(x[,2:100])

# how well did our model do?
sum(x_predict==x_compare)
print(x_predict)
print(x_compare)

# 55! That's reasonably above chance.

nsubs = as.numeric(length(x_all[,1]))

# let's see how the model goes for more than 1 subject. Let's run this on all subjects
pred_success <- array(nsubs)

start_time = Sys.time()

x <- x_all
X <- X_all
ntrials <- ntrials_all

# set up jags and run jags model
data <- list("x","X","ntrials","nsubs") 
params<-c("mu_w","mu_A","mu_theta","mu_a","lambda_w","lambda_A","lambda_theta","lambda_a", "p")
temp_samples <- jags.parallel(data, inits=NULL, params,
                              model.file ="hier_PVL.txt", n.chains=3, 
                              n.iter=3000, n.burnin=1000, n.thin=1, n.cluster=4)

p_post <- temp_samples$BUGSoutput$sims.list$p

post = temp_samples$BUGSoutput$sims.list

for (s in 1:nsubs) {
  
  ntrials = ntrials_all[s] # We don't have predictions for trial 1!
  x_predict <- array(ntrials_all[s])
  
  for (t in 1:(ntrials-1)) { #Trials go from 1 to 99, but correspond to 2-100.
    p_predict <- c(
      MPD(p_post[,s,t,1]),
      MPD(p_post[,s,t,2]),
      MPD(p_post[,s,t,3]),
      MPD(p_post[,s,t,4])
    )
    x_predict[t+1] <- which.max(p_predict) # Increment x-predict by 1 before predicting
    
  }
  print(x_predict[20:40])
  pred_success[s] <- sum(x_predict==x[s, 2:ntrials]) # only comparing with appropriate subject and trials for which we have choices
  print(s)
  
}

end_time = Sys.time()
end_time - start_time

pred_success_adjust <- pred_success/100

save(pred_success_adjust, file = "posterior_predictive_checks.rda")
rm(pred_success_adjust)
load("predictive_checks/posterior_predictive_checks.rda")
avg_pred <- mean(pred_success_adjust) # 0.51. Could be worse.

# plotting code courtesy of Mia
pred_df <- data.frame(pred_success_adjust)
pred_df$sub <- 1:length(pred_success_adjust) # rownames(pred_df) # creating a subject index
pred_df$avg <- mean(pred_df$pred_success)
pred_df$std <- sd(pred_df$pred_success)
pred_df$chance <- .25
ggplot(pred_df, aes(sub, pred_success_adjust)) +
  geom_point() +
  geom_line(aes(y=chance), linetype="dashed", color = "black") +
  geom_ribbon(aes(xmin = -Inf, xmax = Inf, ymin = avg - std, ymax = avg + std), fill = "pink", alpha = 0.6) + 
  geom_line(aes(y=avg))


write.csv(pred_df, "out/pred_success.txt")
