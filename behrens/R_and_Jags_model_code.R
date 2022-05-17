library(tidyverse)
library(rjags)
library(coda)
library(MCMCvis)

# create data
meta = data.frame(probs = c(0.9,0.2, 0.8,0.2,0.8, 0.2, 0.8, 0.1), n = c(120,40,30,40,30,40,40,120))
data = data.frame(reward = unlist(map2(meta$probs, meta$n, function(probs,n){sample(c(0,1),n,T,prob =c(1-probs,probs))})))

# Bayesian learner model code
mod = 'model
{
  R[1] = 0.5
  K ~ dnorm(0.00000E+00, 1/1000)
  V[1] ~ dnorm(0.00000E+00, 1/exp(K))
  for (t in 2:length(reward)) {
    V[t] ~ dnorm(V[t - 1], 1/exp(K))
    R[t] ~ dnorm(R[t - 1], 1/exp(V[t - 1]))  T(0.00000E+00, 
                                               1)
    C[t - 1] <- ifelse(R[t - 1] > 0.5, 1, 0.00000E+00)
    reward[t - 1] ~ dbern(R[t - 1])
  }
}'


# no initial values
inits<-NULL

# what parameters we want to track
params = c("R", "K", "V")

# burn in interval
nb = 30000
# thinning interval
nt = 1
# number of chains
nc = 16


# compile model
jmod = jags.model(file = textConnection(mod), data = data, n.chains = nc, inits = inits, n.adapt = 1000)

# iterate through jmod for the extent of the burn-in
update(jmod, n.iter=nb, by=1)

# draw samples from the posterior for params
post = coda.samples(jmod, params, n.iter = 1000, thin = nt)

#Gelman–Rubin convergence check
gelman.diag(post, multivariate = F)

# diagnostic evaluation of posterior samples
MCMCtrace(post, params = c("R", "K", "V"), pdf=F)


# get posterior means for plot
posterior_means_V = MCMCsummary(post, params = c("V"), digits=2)$mean
posterior_means_R = MCMCsummary(post, params = c("R"), digits=2)$mean


data_for_plot = data.frame(V = posterior_means_V,R =posterior_means_R, probability = unlist(map2(meta$probs, meta$n, function(probs,n){rep(probs, n)})))

# plot estimated probability and real probability
ggplot(data=data_for_plot, mapping = aes(x = 1:nrow(data_to_plot))) + geom_line(mapping = aes(y = R, colour="Estimated probability")) +
  geom_line(mapping = aes(y = probability, colour="Real probability")) +
  scale_color_manual(name = "", values = c("Estimated probability" = "darkblue", "Real probability" = "red")) +
  theme(legend.position = "bottom") +
  xlab("Trial") + 
  ylab("Probability")

# plot estimated volatility
ggplot(data=data_for_plot, mapping = aes(x = 1:nrow(data_to_plot))) + geom_line(mapping = aes(y = V,colour="Volatility")) +  scale_color_manual(name = "", values = c("Volatility" = "darkgreen")) +
  theme(legend.position = "bottom") +
  xlab("Trial") + 
  ylab("Volatility")

