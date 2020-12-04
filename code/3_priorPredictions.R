# McCain et al. 2020
# Unusually large upward shifts in cold-adapted, montane mammals as temperature warms
# 
# Prior predictive checks


############--------------------------------------------------------------------
## Description
##-------------------------------------

# This script runs prior predictive checks to confirm that the prior 
# distributions described in code/models/1_multinom_main.txt are reasonable. 
# Logit-linked regressions can have difficulty converging with overly broad 
# prior distributions, so we chose priors that resulted in plausible range
# extensions of approximately 0-600m beyond the observed range boundary.








############--------------------------------------------------------------------
## Set up
##-------------------------------------

library(truncnorm); library(tidyverse)

inv_logit <- function(x) {1/(1+exp(-x))}

sets <- expand.grid(era=c("H", "C"),
                    mtn=c("FR", "SJ"))
n_sim <- 1e6
priors <- list(alpha=rnorm(n_sim, -10, 1), 
               sigma=rtruncnorm(n_sim, a=0, mean=0, sd=sqrt(1/0.1)),
               beta1=rtruncnorm(n_sim, b=0, mean=-20, sd=sqrt(1/0.5)),
               beta2=rtruncnorm(n_sim, a=0, mean=1, sd=sqrt(1/0.5)))
priors$a <- rnorm(n_sim, priors$alpha, priors$sigma)








############--------------------------------------------------------------------
## Run predictive checks
##-------------------------------------

for(i in 1:nrow(sets)) {
  obs.i <- read_csv(paste0("out/obs/obs_", i, ".csv"))
  
  # assess all combinations of distAway x interpPatchy
  prior_test.df <- obs.i %>% 
    mutate(interpPatchy_sc=scale(interpPatchy)) %>%
    group_by(distAway_sc, interpPatchy_sc) %>% 
    sample_n(1) %>%
    arrange(distAway) %>% ungroup %>% filter(distAway < 1000)
  
  preds <- matrix(NA, nrow(prior_test.df), n_sim)
  
  for(j in 1:n_sim) {
    preds[,j] <- inv_logit(priors$a[j] + 
                             priors$beta1[j]*prior_test.df$distAway_sc + 
                             priors$beta2[j]*prior_test.df$interpPatchy_sc)
  }
  
  
  prior_test.df$pred_mn <- rowMeans(preds)
  prior_test.df$pred_lo <- HDInterval::hdi(t(preds))[1,]
  prior_test.df$pred_hi <- HDInterval::hdi(t(preds))[2,]
  
  pdf(paste0("figs/prior_checks/priors_", sets$mtn[i], "_", sets$era[i], ".pdf"), 
      width=10, height=7)
  par(mfrow=c(2,3))
  plot(density(priors$alpha), xlim=c(-20, 10), main=expression(Prior:~alpha))
  plot(density(priors$beta1), xlim=c(-30, 0), main=expression(Prior:~beta[1]))
  plot(density(priors$beta2), xlim=c(0, 10), main=expression(Prior:~beta[2]))
  plot(pred_lo ~ distAway, data=prior_test.df, ylim=c(0,1), col=rgb(0,0,0,0.2),
       xlab="Distance from empirical range", ylab="Pr(presence)", 
       main=expression(psi[95~HDI~low]))
  plot(pred_mn ~ distAway, data=prior_test.df, ylim=c(0,1), col=rgb(0,0,0,0.2),
       xlab="Distance from empirical range", ylab="Pr(presence)", 
       main=expression(bar(psi)))
  plot(pred_hi ~ distAway, data=prior_test.df, ylim=c(0,1), col=rgb(0,0,0,0.2),
       xlab="Distance from empirical range", ylab="Pr(presence)", 
       main=expression(psi[95~HDI~high]))
  dev.off()
}



