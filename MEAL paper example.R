### *******************************************************************************
##
## Risk Estimation via "Maximum Entropy Acceptable Likelihood" (MEAL) method
## From Faynzilberg (1997), with many adaptations
## by Max Mauerman
##
## goal is to estimate a probability *distribution* of severity for each year
##
## applied to toy data for paper 
##
## this version calculates distributions year-by-year. in general, this should
## be equivalent to the old version, but gets around the need to impute missing data
##
### ********************************************************************************

library(gtools)
library(CVXR)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(dplyr)
library(see)
library(readxl)
library(RColorBrewer)
setwd("/Users/mmauerman/Documents/IRI - Satwin/")

### generate data

prior_data <- c(1,NA,5,2,4)
  
main_data <- cbind(c(1,2,0,5,4),c(1,2,0,5,3),c(1,2,1,0,4))

## user-defined params: data source, number of years (n), number of prob. breaks (l)
## number of error breaks (f), ML threshold (alpha), accuracy / precision tradeoff (beta)
## weight on prior (gamma)

prior_data <- data.frame(prior_data)
data <- data.frame(main_data) 
n <- 5
l <- 5  # higher = greater resolution but more computation needed
f <- 5 # higher = tighter errors but more computation needed
alpha <- 10^-10 # lower = greater weight on ML
beta = 20 # higher = more precision, less accuracy
gamma = 5 # higher = greater weight on prior

### define state-space

# generate S[l] state-space for main model

s <- seq(from = 0, to = 1, by = 1/l)
m <- length(s)

# generate E(f) state-space for error

e <- seq(from = -1, to = 1, by = 1/f)
k <- length(e)

### prepare data

## clean

year <- c(2000:2004)
rownames(data) <- year
rownames(prior_data) <- year

## deal with missing values

# drop years with only missing 

data <- data[rowSums(is.na(data)) != ncol(data), ]

## subset to most recent n years

data <- data[c((nrow(data)-n+1) : nrow(data)),]

# drop vars with only missing or one value

data <- data %>%
  select_if(~sum(!is.na(.)) > 1)

# for ranking purposes, treat non-mentioned years as definitely not bad

prior_data[is.na(prior_data)] <- -99

## convert to quantiles

data <- apply(data,2,function(x) rank(x,ties.method="max",na.last="keep"))
data <- apply(data,2,function(x) x / max(x,na.rm=T))

prior_data <- apply(prior_data,2,function(x) rank(x,ties.method="max",na.last="keep"))
prior_data <- apply(prior_data,2,function(x) x / max(x,na.rm=T))

# drop non-mentioned FBYs for subsequent calculations

prior_data <- apply(prior_data,2,function(x) ifelse(x == min(x),NA,x))

# subset to most recent n years

prior_data <- data.frame(prior_data[c((nrow(prior_data)-n+1) : nrow(prior_data)),])

### error prior

w <- Variable(k, nonneg = TRUE)
e2 <- 10*e^2
obj_errprior <- Maximize(sum(entr(w)))
addup_errprior <- (sum(w) == 1)
consistency_errprior <- (t(10*e) %*% w  == 0)
consistency_errprior_sd <- (t(e2) %*% w  == 1)
prob_errprior <- Problem(obj_errprior,list(addup_errprior,consistency_errprior,
                                           consistency_errprior_sd))
result_errprior <- solve(prob_errprior,solver="SCS")
prior_error <- result_errprior$getValue(w)

### calculate most likely observed states for each predictor ("freq") using MaxEnt 

freq_indiv_eachyear <- list()

for (j in 1:ncol(data)) {
  
  pmf_indiv <- data[,j]
  
  for (y in c(1:n)) {
    
    if (!is.na(pmf_indiv[y])) {
      
      q <- Variable(m, nonneg = TRUE)
      error <- Variable(k, nonneg = TRUE)
      obj <- Maximize(sum(entr(q)) + beta*sum(entr(error)) + beta*sum(error*log(prior_error)))
      addup <- (sum(q) == 1)
      addup_error <- (sum(error) == 1)
      consistency <- (t(s) %*% q + t(e) %*% error == pmf_indiv[y])
      prob_maxent <- Problem(obj,list(addup,addup_error,consistency))
      result_maxent <- solve(prob_maxent,solver = "SCS")
      est_states_freq_indiv_part <- result_maxent$getValue(q)
      
    }
    
    else {
      
      est_states_freq_indiv_part <- rep(NA,m)
      
    }
    
    if (y == 1) {
      
      est_states_freq_indiv <- est_states_freq_indiv_part
      
    }
    
    else {
      
      est_states_freq_indiv <- cbind(est_states_freq_indiv,est_states_freq_indiv_part)  
      
    }
    
  }  
  
  freq_indiv_eachyear[[j]] <- est_states_freq_indiv
}

freq <- apply(simplify2array(freq_indiv_eachyear), 1:2, mean,na.rm=T)

### calculate overall observed means ("pmf")

pmf <- rowMeans(data,na.rm=T)

### prior distribution

prior_indiv_eachyear <- list()

for (j in 1:ncol(prior_data)) {
  
  prior_indiv <- prior_data[,j]
  
  for (y in c(1:n)) {
    
    
    if (!is.na(prior_indiv[y])) {
      
      p <- Variable(m, nonneg = TRUE)
      error <- Variable(k, nonneg = TRUE)
      obj_prior <- Maximize(sum(entr(p)) + beta*sum(entr(error)) + beta*sum(error*log(prior_error)))
      addup_prior <- (sum(p) == 1)
      addup_error <- (sum(error) == 1)
      consistency_prior <- (t(s) %*% p + t(e) %*% error == prior_indiv[y])
      prob_prior <- Problem(obj_prior,list(addup_prior,addup_error,consistency_prior))
      result_prior <- solve(prob_prior,solver="SCS")
      
      est_states_prior_indiv <- result_prior$getValue(p)
      
    }
    
    else {
      
      est_states_prior_indiv <- rep(1/m,m)
      
    }
    
    if (y == 1) {
      
      prior_dist_indiv <- est_states_prior_indiv
      
    }
    
    else {
      
      prior_dist_indiv <- cbind(prior_dist_indiv,est_states_prior_indiv)
      
    }
    
  }
  
  prior_indiv_eachyear[[j]] <- prior_dist_indiv
  
}

prior <- apply(simplify2array(prior_indiv_eachyear), 1:2, mean,na.rm=T)

### solve main model

for (y in c(1:n)) {
  
  ### define optimization conditions
  
  q <- Variable(m, nonneg = TRUE)
  error <- Variable(k, nonneg = TRUE)
  
  ## entropy
  
  obj <- Maximize(sum(entr(q)))
  obj_withprior <- Maximize(sum(entr(q)) + gamma*sum(q*log(prior[,y])) +
                              beta*sum(entr(error)) + beta*sum(error*log(prior_error)))
  
  ## adding-up constraints
  
  addup <- (sum(q) == 1)
  addup_error <- (sum(error) == 1)
  
  ## aggregate data consistency
  
  consistency <- (t(s) %*% q + t(e) %*% error == pmf[y])
  consistency_noerror <- (t(s) %*% q == pmf[y])
  
  ## maximum likelihood function 
  
  nonmiss <- sum(!is.na(data[y,]))
  ml <- (sum(freq[,y] * (log(freq[,y]) - log(q))) <= (1/nonmiss)*log(1/alpha))
  
  ### solve 
  
  prob <- Problem(obj_withprior,list(addup,addup_error,consistency,ml))
  result <- solve(prob,solver="SCS")
  est_states_indiv <- result$getValue(q)
  est_error_indiv <- result$getValue(error)
  
  ## solve without MLE or prior for comparison 
  
  prob_maxent <- Problem(obj,list(addup,consistency_noerror))
  result_maxent <- solve(prob_maxent,solver="SCS")
  est_states_maxent_indiv <- result_maxent$getValue(q)
  
  if (y == 1) {
    
    est_states <- est_states_indiv
    est_states_maxent <- est_states_maxent_indiv
    est_error <- est_error_indiv
    
  }
  
  else {
    
    est_states <- cbind(est_states,est_states_indiv)
    est_states_maxent <- cbind(est_states_maxent,est_states_maxent_indiv)
    est_error <- cbind(est_error,est_error_indiv)
    
  }
  
}

### aggregate and plot

## calculate estimated prob dist for each year 

# main model 

est_probs <- data.frame(est_states,s)
colnames(est_probs) <- c(rownames(data),"s")
est_probs_long <- gather(est_probs,year,prob,-s)

# maxent only (for comparison)

est_probs_maxent <- data.frame(est_states_maxent,s)
colnames(est_probs_maxent) <- c(rownames(data),"s")
est_probs_long_maxent <- gather(est_probs_maxent,year,prob,-s)

# MLE only (for comparison)

est_probs_mle <- data.frame(freq,s)
colnames(est_probs_mle) <- c(rownames(data),"s")
est_probs_long_mle <- gather(est_probs_mle,year,prob,-s)

# prior only

est_probs_prior <- data.frame(prior,s)
colnames(est_probs_prior) <- c(rownames(data),"s")
est_probs_long_prior <- gather(est_probs_prior,year,prob,-s)

# error term

est_probs_error <- data.frame(est_error,e)
colnames(est_probs_error) <- c(rownames(data),"e")
est_probs_long_error <- gather(est_probs_error,year,prob,-e)

# uniform dist (for reference)

est_states_unif <- rep(1/length(s),length(s))
est_states_unif <- matrix(rep(as.numeric(est_states_unif),each=n),nrow=length(est_states_unif))

est_probs_unif <- data.frame(est_states_unif,s)
colnames(est_probs_unif) <- c(rownames(data),"s")
est_probs_long_unif <- gather(est_probs_unif,year,prob,-s)

### Plots

# main

# est_probs_long_nozero <- est_probs_long[est_probs_long$prob != 0,]

plot <- ggplot(est_probs_long, aes(x=factor(year),group=factor(s),y=prob,fill=s)) +
  geom_col(position=position_dodge()) + ylim(0,1) +
  scale_fill_continuous(type="viridis") +
  theme_minimal() + theme(legend.position = "none") +
  xlab("Year") + ylab("Probability") + ggtitle("MEAL estimate") + coord_flip()

# maxent only

# est_probs_long_maxent_nozero <- est_probs_long_maxent[est_probs_long_maxent$prob != 0,]

plot_maxent <- ggplot(est_probs_long_maxent, aes(x=factor(year),group=factor(s),y=prob,fill=s)) +
  geom_col(position=position_dodge()) + ylim(0,1) +
  scale_fill_continuous(type="viridis") +
  theme_minimal() + theme(legend.position = "none") +
  xlab("Year") + ylab("Probability") + ggtitle("MaxEnt only estimate") + coord_flip()

# MLE only 

# est_probs_long_mle_nozero <- est_probs_long_mle[est_probs_long_mle$prob != 0,]

plot_mle <- ggplot(est_probs_long_mle, aes(x=factor(year),group=factor(s),y=prob,fill=s)) +
  geom_col(position=position_dodge()) + ylim(0,1) +
  scale_fill_continuous(type="viridis") +
  theme_minimal() + theme(legend.position = "none") +
  xlab("Year") + ylab("Probability") + ggtitle("MLE only estimate") + coord_flip()

# prior only 

# est_probs_long_prior_nozero <- est_probs_long_prior[est_probs_long_prior$prob != 0,]

plot_prior <- ggplot(est_probs_long_prior, aes(x=factor(year),group=factor(s),y=prob,fill=s)) +
  geom_col(position=position_dodge()) + ylim(0,1) +
  scale_fill_continuous(type="viridis") +
  theme_minimal() + theme(legend.position = "none") +
  xlab("Year") + ylab("Probability") + ggtitle("prior only estimate") + coord_flip()

# all together

grid.arrange(plot_maxent,plot_mle,plot_prior,plot,nrow=1,ncol=4)

## diagnostics

# expected error

for (y in c(1:n)) {
  
  est_expval_part <- est_states[,y] * s
  
  if (y == 1) {
    
    est_expval <- est_expval_part
    
  }
  
  else {
    
    est_expval <- cbind(est_expval,est_expval_part)
    
  }
  
}

est_expval <- colSums(est_expval)
est_expval <- data.frame(cbind(est_expval,rownames(data)))
colnames(est_expval) <- c("expval","year")

e_error <- est_probs_long_error %>%
  group_by(year) %>%
  summarise(e_error = sum((prob*e))) %>%
  mutate(e_error_overall = mean(abs(e_error)))

e_error <- left_join(e_error,est_expval,by="year")
e_error$expval <- as.numeric(e_error$expval)

# Efficiency (information content)

for (y in c(1:n)) {
  
  efficiency_part <- -sum(est_states[,y] * log(est_states[,y])) / -sum(est_states_unif * log(est_states_unif))
  
  if (y == 1) {
    
    efficiency_byyear <- efficiency_part
    
  }
  
  else {
    
    efficiency_byyear <- c(efficiency_byyear,efficiency_part)
    
  }
  
}

efficiency <- mean(efficiency_byyear,na.rm=T)

# cross-entropy

for (y in c(1:n)) {
  
  crossentr_part <- -sum((est_states[,y] * log(prior[,y]))) / -sum(prior[,y] * log(prior[,y]))
  
  if (y == 1) {
    
    crossentr_byyear <- crossentr_part
    
  }
  
  else {
    
    crossentr_byyear <- c(crossentr_byyear,crossentr_part)
    
  }
  
}

crossentr_byyear <- data.frame(cbind(crossentr_byyear,rownames(data)))
colnames(crossentr_byyear) <- c("crossentr","year")
crossentr_byyear$crossentr <- as.numeric(crossentr_byyear$crossentr)
crossentr_byyear$crossentr <- format(crossentr_byyear$crossentr,digits=2,nsmall=2)

crossentr_lab <- paste(crossentr_byyear$year,'\n cross-entropy = ',crossentr_byyear$crossentr)
names(crossentr_lab) <- as.character(crossentr_byyear$year)

## fancier plots

data_rug <- as.data.frame(t(data))
data_rug <- data_rug %>%
  gather(key="year",value="value")

prior_rug <- as.data.frame(t(prior_data))
prior_rug <- prior_rug %>%
  gather(key="year",value="value")

pmf_rug <- as.data.frame(t(pmf))
pmf_rug <- pmf_rug %>%
  gather(key="year",value="value")

newplot <- ggplot(est_probs_long,aes(x=s,y=prob)) +
  geom_col(position=position_nudge(x=-0.0),fill="red") +
  geom_rug(data=data_rug,aes(x=value),inherit.aes = F) +
  geom_rug(data=prior_rug,aes(x=value),inherit.aes = F,color="blue",size=1) +
  geom_rug(data=pmf_rug,aes(x=value),inherit.aes = F,color="green",size=1.5) +
  geom_col(data=est_probs_long_unif,position=position_nudge(x=-0.0),aes(x=s,y=prob),alpha = 0.3) +
  ylim(0,1) +
  facet_wrap(~year,labeller = labeller(year=crossentr_lab)) +
  theme_minimal() +
  xlab("Severity state") + ylab("Probability of state") + 
  ggtitle("Estimated severity by year",subtitle="X axis = return period. Grey bars are uniform dist. for reference. Observations in margins; prior data in blue")
newplot

est_probs_long_comparison <- left_join(est_probs_long,est_probs_long_unif,by=c("s" = "s","year" = "year"))
est_probs_long_comparison$diff <- est_probs_long_comparison$prob.x - est_probs_long_comparison$prob.y

newplot_comparison <- ggplot(est_probs_long_comparison,aes(x=s,y=diff)) +
  geom_col(position=position_nudge(x=-0.0),fill="red") +
  geom_rug(data=data_rug,aes(x=value),inherit.aes = F) +
  geom_rug(data=prior_rug,aes(x=value),inherit.aes = F,color="blue",size=1) +
  geom_rug(data=pmf_rug,aes(x=value),inherit.aes = F,color="green",size=1.5) +
  ylim(-0.75,0.75) +
  facet_wrap(~year,labeller = labeller(year=crossentr_lab)) +
  theme_minimal() +
  xlab("Severity state") + ylab("Relative probability of state") + 
  ggtitle("Estimated severity by year",subtitle="X axis = return period. Observations in margins; prior data in blue")
newplot_comparison

est_probs_long_exceed <- est_probs_long %>%
  group_by(year) %>%
  mutate(cumprob = sapply(s,function(x) sum(prob[s >= x]) ))

est_probs_long_unif_exceed <- est_probs_long_unif %>%
  group_by(year) %>%
  mutate(cumprob = sapply(s,function(x) sum(prob[s >= x]) ))

est_probs_long_exceed$cumprob[est_probs_long_exceed$cumprob < 0] <- 0

newplot_exceed <- ggplot(est_probs_long_exceed,aes(x=s,y=cumprob)) +
  geom_col(position=position_nudge(x=0.0),fill="red") +
  geom_rug(data=data_rug,aes(x=value),position=position_nudge(x=-0.0),inherit.aes = F) +
  geom_rug(data=prior_rug,aes(x=value),position=position_nudge(x=-0.0),inherit.aes = F,color="blue",size=1) +
  geom_rug(data=pmf_rug,aes(x=value),position=position_nudge(x=-0.0),inherit.aes = F,color="green",size=1.5) +
  geom_col(data=est_probs_long_unif_exceed,position=position_nudge(x=0.0),aes(x=s,y=cumprob),alpha = 0.3) +
  ylim(0,1.1) +
  facet_wrap(~year,labeller = labeller(year=crossentr_lab)) +
  theme_minimal() +
  xlab("Severity state") + ylab("Probability of exceedence") + 
  ggtitle("Estimated severity by year",subtitle="X axis = return period. Grey bars are uniform dist. for reference. Observations in margins; prior data in blue")
newplot_exceed

est_probs_long_nonexceed <- est_probs_long %>%
  group_by(year) %>%
  mutate(cumprob = sapply(s,function(x) sum(prob[s <= x]) ))

est_probs_long_unif_nonexceed <- est_probs_long_unif %>%
  group_by(year) %>%
  mutate(cumprob = sapply(s,function(x) sum(prob[s <= x]) ))

est_probs_long_nonexceed$cumprob[est_probs_long_nonexceed$cumprob < 0] <- 0

newplot_nonexceed <- ggplot(est_probs_long_nonexceed,aes(x=s,y=cumprob)) +
  geom_col(position=position_nudge(x=0.0),fill="red") +
  geom_rug(data=data_rug,aes(x=value),position=position_nudge(x=-0.0),inherit.aes = F) +
  geom_rug(data=prior_rug,aes(x=value),position=position_nudge(x=-0.0),inherit.aes = F,color="blue",size=1) +
  geom_rug(data=pmf_rug,aes(x=value),position=position_nudge(x=-0.0),inherit.aes = F,color="green",size=1.5) +
  geom_col(data=est_probs_long_unif_nonexceed,position=position_nudge(x=0.0),aes(x=s,y=cumprob),alpha = 0.3) +
  ylim(0,1.1) +
  facet_wrap(~year,labeller = labeller(year=crossentr_lab)) +
  theme_minimal() +
  xlab("Severity state") + ylab("Probability of non-exceedence") + 
  ggtitle("Estimated severity by year",subtitle="X axis = return period. Grey bars are uniform dist. for reference. Observations in margins; prior data in blue")
newplot_nonexceed

est_probs_long_error_labs <- est_probs_long_error[est_probs_long_error$prob != 0,]
est_probs_long_error_labs <- left_join(est_probs_long_error,e_error,by=c("year" = "year"))

e_error_formatted <- format(max(e_error$e_error_overall),digits=2,nsmall=2)

plot_error <- ggplot(est_probs_long_error_labs, aes(x=factor(year),group=factor(e),y=prob,fill=e)) +
  geom_col(position=position_dodge()) + ylim(0,0.3) +
  scale_fill_continuous(type="gradient") +
  theme_minimal() + theme(legend.position = "none") +
  geom_label(aes(label=format(e_error,digits=2,nsmall=2),y=0.3,x=factor(year)),inherit.aes = FALSE) +
  xlab("Year") + ylab("Probability of error state") + labs(title = "Error term",
                                                           subtitle = "Extreme X values = larger errors. Expected error in boxes",
                                                           caption=paste("Model efficiency =",sprintf("%.1f %%", 100*(1-efficiency)),", Average error = ",e_error_formatted))
plot_error





