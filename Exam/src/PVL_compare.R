install.packages("pacman")
pacman::p_load(pacman,R2jags, parallel, polspline, tidyverse)
setwd("/work/Exam")
# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#load control data
choices <- read.csv("raw_data/choice.csv",header=TRUE)
outcomes <- read.csv("raw_data/outcome.csv",header=TRUE)
demographics <- read.csv("preprocessed/demographics.csv",header=TRUE)
demographics = subset(demographics, select = c("ID", "Gender"))
## SANITY CHECK - ID colnames appear misspecified by 155.
sanity_1 <- read.csv("raw_data/Gamtest_PreprocessedOSF.csv")
sanity_2 <- read.csv("raw_data/choice.csv")
sanity_names_1 <- colnames(sanity_2[,2:141])
sanity_names_2 <- c()

for (name in c(sanity_1$ID - 155)){
  string = paste('X', as.character(name), sep = "")
  sanity_names_2 = c(sanity_names_2, string)
}
print(sum(sanity_names_1 == sanity_names_2)) # Perfect overlap.
## SANITY CHECK OVER
## THEREFORE:
demographics$ID = demographics$ID - 155 # See sanity check below for why this is apparantly necessary.


# Get corresponding ID's
female_IDs = c()
male_IDs = c()
demographics_male = filter(demographics, Gender == "Male")
demographics_female = filter(demographics, Gender == "Female")

for (name in demographics_male$ID){
  string = paste('X', as.character(name), sep = "")
  male_IDs = c(male_IDs, string)

}

for (name in demographics_female$ID){
  string = paste('X', as.character(name), sep = "")
  female_IDs = c(female_IDs, string)
  
}
# Saving and SCALING!!! variables
choices_men <- subset(choices, select = male_IDs)
choices_women <- subset(choices, select = female_IDs)
outcomes_men <- subset(outcomes, select = male_IDs)
outcomes_women <- subset(outcomes, select = female_IDs)

write.csv(choices_men, "preprocessed/choices_men.csv")
write.csv(choices_women, "preprocessed/choices_women.csv")
write.csv(outcomes_men, "preprocessed/outcomes_men.csv")
write.csv(outcomes_women, "preprocessed/outcomes_women.csv")





## Run from here if above chunk has been run (remember to load packages and define functions first!)
set.seed(1337)
choices_men <- read.csv("preprocessed/choices_men.csv")[,2:19]
choices_women <- read.csv("preprocessed/choices_women.csv")[,2:21]
outcomes_men <- read.csv("preprocessed/outcomes_men.csv")[,2:19]/100
outcomes_women <- read.csv("preprocessed/outcomes_women.csv")[,2:21]/100


x_men = as.data.frame(t(choices_men))
X_men = as.data.frame(t(outcomes_men))

x_women = as.data.frame(t(choices_women))
X_women = as.data.frame(t(outcomes_women))

nsubs_men = as.numeric(length(choices_men[1,]))
nsubs_women = as.numeric(length(choices_women[1,]))

ntrials <- 100


# set up jags and run jags model on one subject
data <- list("x_men","X_men","nsubs_men", "ntrials",
             "x_women","X_women", "nsubs_women") 
params<-c("alpha_a","alpha_A","alpha_theta","alpha_w", "mu_a", "mu_A", "mu_theta", "mu_w")

start_time = Sys.time()
samples <- jags.parallel(data, inits=NULL, params,
                model.file ="compare/PVL_compare.txt",
                n.chains=3, n.iter=30000, n.burnin=1000, n.thin=2, n.cluster=4)

end_time = Sys.time()
end_time - start_time

# Save results and environment!
save.image(file = "out/results_environment.RData")


# Traceplots and overall output:
main_analysis = samples$BUGSoutput$summary
write.table(samples$BUGSoutput$summary, "out/main_analysis.txt",sep="\t",row.names=TRUE)

traceplot(samples)

BF_effect <- NULL
BF_null <- NULL
MPDs <- NULL
#### - alpha_a
# savage dickey plot - alpha_a
plot(density(rnorm(12000,0.2,1)),,ylim=c(0,3),main="alpha_a")

lines(density(samples$BUGSoutput$sims.list$alpha_a),col="red")

# Maximum posterior density
MPDs$alpha_a <- MPD(samples$BUGSoutput$sims.list$alpha_a)

#Logsplines for Bayes Factor calculations
fit.posterior <- logspline(samples$BUGSoutput$sims.list$alpha_a)
effect.posterior <- dlogspline(0.2, fit.posterior)
effect.prior <- dnorm(0.2,0.2,1)

# Logsplines for zero calculations
null.posterior <- dlogspline(0, fit.posterior)
null.prior     <- dnorm(0,0,1)                   

BF_effect$alpha_a <- effect.posterior/effect.prior # Posterior is approx. 0.33
# Aka, one might say that we used to believe 3 times as much in 0.2 as we do now.

BF_null$alpha_a<- null.posterior/null.prior


############## alpha_w #############
# savage dickey plot - alpha_w
plot(density(rnorm(12000,0,1/sqrt(1))),ylim=c(0,3),main="alpha_w")

lines(density(samples$BUGSoutput$sims.list$alpha_w),col="red")

# Maximum posterior density
MPDs$alpha_w <- MPD(samples$BUGSoutput$sims.list$alpha_w)

#Logsplines for Bayes Factor calculations
fit.posterior <- logspline(samples$BUGSoutput$sims.list$alpha_w)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))

BF_null$alpha_w<- null.posterior/null.prior
BF_effect$alpha_w <- null.prior/null.posterior



############## alpha_A #############
# savage dickey plot - alpha_A
plot(density(rnorm(12000,0,1/sqrt(1))),ylim=c(0,4),main="alpha_A")

lines(density(samples$BUGSoutput$sims.list$alpha_A),col="red")

# Maximum posterior density
MPDs$alpha_A <- MPD(samples$BUGSoutput$sims.list$alpha_A)

#Logsplines for Bayes Factor calculations
fit.posterior <- logspline(samples$BUGSoutput$sims.list$alpha_A)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))

BF_null$alpha_A<- null.posterior/null.prior
BF_effect$alpha_A <- null.prior/null.posterior


############## alpha_theta #############
# savage dickey plot - alpha_theta
plot(density(rnorm(12000,0,1/sqrt(1))),ylim=c(0,3),main="alpha_theta")

lines(density(samples$BUGSoutput$sims.list$alpha_theta),col="red")

# Maximum posterior density
MPDs$alpha_theta <- MPD(samples$BUGSoutput$sims.list$alpha_theta)

#Logsplines for Bayes Factor calculations
fit.posterior <- logspline(samples$BUGSoutput$sims.list$alpha_theta)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))

BF_null$alpha_theta<- null.posterior/null.prior
BF_effect$alpha_theta <- null.prior/null.posterior


############## mu_a #############
# savage dickey plot - mu_a
plot(density(rnorm(12000,0,1/sqrt(1))),ylim=c(0,4),main="mu_a")

lines(density(samples$BUGSoutput$sims.list$mu_a),col="red")

# Maximum posterior density
MPDs$mu_a <- MPD(samples$BUGSoutput$sims.list$mu_a)


############## mu_w #############
# savage dickey plot - mu_w
plot(density(rnorm(12000,0,1/sqrt(1))),ylim=c(0,3),main="mu_w")

lines(density(samples$BUGSoutput$sims.list$mu_w),col="red")

# Maximum posterior density
MPDs$mu_w <- MPD(samples$BUGSoutput$sims.list$mu_w)


############## mu_theta #############
# savage dickey plot - mu_theta
plot(density(rnorm(12000,0,1/sqrt(1))),ylim=c(0,7),main="mu_theta")

lines(density(samples$BUGSoutput$sims.list$mu_theta),col="red")

# Maximum posterior density
MPDs$mu_theta <- MPD(samples$BUGSoutput$sims.list$mu_theta)


############## mu_A #############
# savage dickey plot - mu_A
plot(density(rnorm(12000,0,1/sqrt(1))),ylim=c(0,3),main="mu_A")

lines(density(samples$BUGSoutput$sims.list$mu_A),col="red")

# Maximum posterior density
MPDs$mu_A <- MPD(samples$BUGSoutput$sims.list$mu_A)


############### Look at cumulative balance to justify learning ##############
set.seed(1337)
# calculate cumulative value
xcum_men <- as.data.frame(array(0,c(nsubs_men,100)))

for (i in 1:nsubs_men) {
  xcum_men[i,] <- cumsum(as.numeric(X_men[i,])) 
  print(i)
}

xcum_women <- as.data.frame(array(0,c(nsubs_women,100)))
for (i in 1:nsubs_women) {
  xcum_women[i,] <- cumsum(as.numeric(X_women[i,])) 
}
# Scaling variables
xcum_men = xcum_men
xcum_women = xcum_women

plot(colMeans(xcum_men),ylim=c(-15,8),type='l',lwd=2,  
     xlab = "Trial", ylab = "Cumulative Balance (hundreds of DKK)")

lines(colMeans(xcum_women),col="red",lwd=2)

# Indicates that a useful prior could be a mean around 0 with a standard deviation of 2.

############### compare final means ##############

#sd of 5 = 1/sqrt(1/25)

# set up jags and run jags model on one subject
data <- list("xcum_men","nsubs_men",
             "xcum_women","nsubs_women") 
params<-c("mu","alpha")
temp_samples <- jags.parallel(data, inits=NULL, params,
                model.file ="compare/gender_balance_compare.txt",
                n.chains=3, n.iter=15000, n.burnin=1000, n.thin=3, n.cluster=4)


# Traceplots and overall output:
outcome_analysis = temp_samples$BUGSoutput$summary
traceplot(temp_samples$BUGSoutput)



write.table(outcome_analysis, "out/sanity_check_analysis.txt",sep="\t",row.names=TRUE)



############## alpha #############
# savage dickey plot - alpha
plot(density(rnorm(12000,10,1/sqrt(0.04))), ylim=c(0,0.4),main="alpha")

lines(density(temp_samples$BUGSoutput$sims.list$alpha),col="red")

# Maximum posterior density
MPDs$alpha <- MPD(temp_samples$BUGSoutput$sims.list$alpha)

#Logsplines for Bayes Factor calculations
fit.posterior <- logspline(temp_samples$BUGSoutput$sims.list$alpha)
null.posterior <- dlogspline(0, fit.posterior)
null.prior     <- dnorm(0,10,(1/sqrt(0.04)))
effect.posterior <- dlogspline(10, fit.posterior)
effect.prior <- dnorm(10,mean = 10,sd = 1/sqrt(0.04))

BF_null$alpha<- null.posterior/null.prior
BF_effect$alpha <- effect.posterior/effect.prior




############## mu #############
# savage dickey plot - mu
plot(density(rnorm(12000,0,1)),ylim=c(0,3),main="mu")

lines(density(temp_samples$BUGSoutput$sims.list$mu),col="red")

# Maximum posterior density
MPDs$mu <- MPD(temp_samples$BUGSoutput$sims.list$mu)

write.table(BF_null, "out/BF_null.txt",sep="\t",row.names=TRUE)
write.table(BF_effect, "out/BF_effect.txt",sep="\t",row.names=TRUE)
write.table(MPDs, "out/MPDs.txt",sep="\t",row.names=TRUE)

save.image(file = "FullResults.RData")
