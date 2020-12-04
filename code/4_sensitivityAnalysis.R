# McCain et al. 2020
# Unusually large upward shifts in cold-adapted, montane mammals as temperature warms
# 
# Sensitivity Analysis


############--------------------------------------------------------------------
## Description
##-------------------------------------

# This script simulates communities, then simulates presence-only samples 
# mimicking the pattern of observations in the empirical datasets. It varies
# the size of the elevational bin (50m, 100m, 200m), as well as the sampling
# effort, which is implemented as a proportion of the empirical sampling effort
# on each mountain in each era. See Appendix S2: Sensitivity Analysis for full 
# descriptions of the methods.
#
# The community at each elevation is simulated using a Ricker model for each
# species. The mean intrinsic growth rate for each species is determined by its 
# relationship with elevation, such that for species j at elevational i:
#
#     r.mn[i,j] = b0[j] + b1[j]*el[i] + b2[j]*el[i]^2
#
# And the growth rate in each year k is stochastic:
#
#     r[i,j,k] ~ Norm(r.mn[i,j]*(1 - N[i,j,k-1,]/K[j]), r.sd.[j]) 
#
# with species-specific standard deviation and carrying capacity. 
#
# The species-specific b parameters are distributed normally about global 
# hyperparameters beta, which are set to create a mid-elevation peak in 
# species richness as is observed on both mountains in this study.








############--------------------------------------------------------------------
## Set up
##-------------------------------------

library(tidyverse); library(rjags); library(ggmcmc); library(dclone)
source("code/00_fn.R")

#--- load empirical data
sp.i <- read_csv("data/Mamm_Summary_Data.csv") 
pDet.full <- read_csv("data/prDet_processed.csv")

sets <- expand.grid(b=c(200, 100, 50),
                    mtn=c("FR", "SJ"), 
                    effort=c(0.5, 1, 2),
                    era=c("H", "C"))

#--- set simulation parameters
n.sim <- 50
tmax <- 200  # number of years to simulate
yrs.obs <- (-29:0)+tmax  # years to draw samples from

#--- set hyperparameters
beta <- c(0.2, -0.3, -0.3)  # mean slopes (avg b0, b1, b2 among species)
beta.sd <- c(0.2, 0.5, 0.15)  # sd in slopes (sd of b0, b1, b2 among species)







############--------------------------------------------------------------------
## Determine richness and sampling distribution
##-------------------------------------

for(set in 1:nrow(sets)) {
  
  # parameters for 'set' details
  b <- sets$b[set]  
  mtn <- sets$mtn[set]
  effort <- sets$effort[set]
  era <- sets$era[set]
  
  # subset raw empirical data to era x mtn range
  Y.real <- read_csv(paste0("data/sample_els_", era, ".csv")) %>%
    mutate(elBin=el %/% b * b) 
  if(mtn=="FR") {
    Y.real <- Y.real %>% 
      filter(county %in% c("Boulder", "Larimer", "Big Thompson", "Cedar Park",
                           "Cow Camp", "Denver", "SylvanDale"))
  } else if(mtn=="SJ") {
    Y.real <- Y.real %>%
      filter(county %in% c("Dolores", "La Plata", "Montezuma", "San Juan",
                           "East San Juans", "West San Juans"))
  }
  
  # determine observations per elevational bin
  Y.counts <- Y.real %>% group_by(elBin) %>% summarise(Y=round(n()*effort))
  els <- seq(min(Y.counts$elBin), max(Y.counts$elBin), by=b)
  n.els <- length(els)
  Ytot <- Y.counts$Y[match(els, Y.counts$elBin)] %>% replace_na(0)
  J <- n_distinct(Y.real$sp)
  pDet.set <- pDet.full %>% 
    filter(Species %in% filter(sp.i, Abbrev %in% Y.real$sp)$Species)
  
  
  
  
  ############--------------------------------------------------------------------
  ## Simulate communities, sample, and fit JAGS model
  ##-------------------------------------
  
  gg.Z <- rmse.ls <- vector("list", n.sim)
  
  comm.true <- simulate_communities(J, els, pDet.set, beta, beta.sd, tmax)
  
  #--- example log abundance distribution
  jpeg(paste0("figs/sims/logN_", mtn, "_", era, "_", b, "m.jpg"), 
       width=10, height=5, units="in", res=300)
  {
    par(mfrow=c(1,2)) 
    plot(density(rowSums(log(1+round(comm.true$N[175,,])))), lwd=2,
         main=paste0("Log abundance across gradient:\n",
                     mtn, "-", era, ", ", b, "m bins, year ", 175),
         xlab="log N", ylab="Density")
    legend("topleft", title="A.", legend="", bty="n")
    plot(NA, NA, xlim=c(0, 16), ylim=c(0, 0.3), 
         main=paste0("Log abundance by elevation:\n",
                     mtn, "-", era, ", ", b, "m bins, year ", 175),
         xlab="log N", ylab="Density")
    legend("topleft", title="B.                  ", 
           "Elevation", paste0(els[c(1,n.els/2,n.els)], "m"), bty="n",
           lty=1, lwd=1.5, col=viridis::viridis(n.els+1)[c(1,n.els/2,n.els)])
    for(j in 1:n.els) {
      lines(density(log(round(comm.true$N[175,,j]))), 
            col=viridis::viridis(n.els+1)[j], lwd=1.5)
    }
  }
  dev.off()
  
  
  for(s in 1:n.sim) {
    
    obs <- sample_community(J, els, yrs.obs, comm.true, Ytot)
    
    jags_d <- list(J=J, 
                   n.el=n.els, 
                   el=c(scale(els)),
                   y=obs$Y, 
                   Y=rowSums(obs$Y), 
                   delta_shp=as.matrix(comm.true$pDet.par), 
                   LAMBDA=obs$spAbund,
                   interpPatchy=c(scale(obs$interpPatchy)),
                   distAway=matrix(scale(c(obs$binsAway*b)), ncol=J))
    
    #--- fit model
    pars <- c("Z")
    cl <- makeCluster(4)
    out <- jags.parfit(cl=cl, data=jags_d, params=pars, 
                       model="code/models/1_hb_undersampling.txt",
                       inits=list(Z=matrix(1, n.els, J)),
                       n.chains=3, n.adapt=50000, n.update=20000, 
                       n.iter=20000, n.thin=50)
    stopCluster(cl)
    
    
    ##-- aggregate output
    true.df <- data.frame(bin=rep(1:n.els, times=J),
                          Elevation=rep(els, times=J),
                          spp=rep(1:J, each=n.els),
                          N=c(t(apply(comm.true$N[yrs.obs,,], 2:3, mean))),
                          Z=c(t(apply(comm.true$Z[yrs.obs,,], 2:3, mean)))>0.05,
                          Y=c(obs$Y),
                          Ytot=rep(rowSums(obs$Y), times=J),
                          binsAway=c(obs$binsAway),
                          distAway=c(obs$binsAway)*b,
                          distAway_sc=c(scale(c(obs$binsAway)*b)),
                          interpPatchy=rep(obs$interpPatchy, each=n.els),
                          interp.rng.obs=c(obs$interp.rng),
                          interp.rng.true=c(comm.true$interp.rng.true),
                          spAbund=rep(obs$spAbund, each=n.els),
                          delta=c(t(apply(comm.true$pDet.ar[yrs.obs,,], 
                                          2:3, mean))))
    gg.Z[[s]] <- ggs(out, "Z") %>%
      mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
               str_remove("Z\\[") %>% as.numeric,
             spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
               str_remove("\\]") %>% as.numeric) %>%
      group_by(bin, spp) %>%
      summarise(mn_prPres=mean(value)) %>%
      full_join(true.df, c("bin", "spp")) 
    
    comp.ls[[s]] <- data.frame(spp=1:J, 
                               true.lo=comm.true$rng.med[,1],
                               true.hi=comm.true$rng.med[,2],
                               obs.lo=na_if(obs$rng[1,], 0),
                               obs.hi=na_if(obs$rng[2,], 0)) %>%
      full_join(gg.Z[[s]] %>% filter(mn_prPres >= 0.05) %>%
                  arrange(spp, bin) %>% group_by(spp) %>%
                  summarise(mod.lo=first(bin), 
                            mod.hi=last(bin)),
                by="spp") %>%
      mutate_at(2:7, ~.*b+min(els)) %>% 
      mutate(sim=s, b=b, mtn=mtn, effort=effort, era=era)
    cat("--------------\n", 
        "Finished", s, "of", n.sim, "in set", set, "/", nrow(sets),
        "\n--------------\n")
  }
  
  set_pars <- paste(paste0(names(sets), c(as.character(b), 
                                          as.character(mtn), 
                                          as.character(effort*100),
                                          as.character(era))), 
                    collapse="_")
  
  write_csv(do.call('rbind', comp.ls), 
            paste0("out/sim/comp_", set_pars, ".csv"))
  write_csv(do.call('rbind', gg.NZ),
            paste0("out/sim/ggNZ_", set_pars, ".csv"))
}








############--------------------------------------------------------------------
## Calculate RMSE and accuracy
##-------------------------------------

comp <- map_dfr(dir("out/sim", "comp", full.names=T), read_csv)


#--- RMSE
RMSE.df <- comp %>% filter(!is.na(obs.lo)) %>%
  group_by(b, mtn, effort, era, sim) %>%
  summarise(obs.lo.rmse=sqrt(mean((obs.lo-true.lo)^2, na.rm=T)),
            mod.lo.rmse=sqrt(mean((mod.lo-true.lo)^2, na.rm=T)),
            obs.hi.rmse=sqrt(mean((obs.hi-true.hi)^2, na.rm=T)),
            mod.hi.rmse=sqrt(mean((mod.hi-true.hi)^2, na.rm=T))) %>%
  mutate(diff.lo=(mod.lo.rmse-obs.lo.rmse)/obs.lo.rmse,
         diff.hi=(mod.hi.rmse-obs.hi.rmse)/obs.hi.rmse) %>%
  select(1:5,10:11) %>% ungroup %>%
  pivot_longer(6:7, names_to="boundary", values_to="diff") %>%
  mutate(boundary=factor(str_sub(boundary, -2L, -1L), levels=c("lo", "hi"),
                         labels=c("Lower boundary", "Upper boundary")), 
         b=factor(b, levels=c("25", "50", "100", "200")),
         effort=factor(effort, levels=c(0.5, 1, 2), 
                       labels=c("50%", "100%", "200%")),
         era=factor(era, levels=c("H", "C"), 
                    labels=c("Historical", "Contemporary")))
RMSE.df %>%
  group_by(b, boundary, mtn, effort, era) %>%
  summarise(mnDiff=mean(diff, na.rm=T), 
            seDiff=sd(diff, na.rm=T)/sqrt(max(sim))) %>%
  ggplot(aes(x=effort, y=mnDiff, ymin=mnDiff-2*seDiff, ymax=mnDiff+2*seDiff,
             colour=b, shape=mtn)) +
  geom_hline(yintercept=0, linetype=2) + 
  geom_point(position=position_dodge(width=0.25), size=1.5) + 
  geom_linerange(position=position_dodge(width=0.25)) + 
  facet_grid(era~boundary) + 
  scale_shape_manual("Mountain\nRange", values=c(1,5)) +
  scale_y_continuous(labels=scales::percent) +
  scale_colour_brewer("Bin size", type="qual", palette=2) +
  labs(x="Sampling effort relative to empirical", 
       y="Change in RMSE relative to empirical (mean ± 2 SE)") + 
  theme_bw()
ggsave("figs/sims/RMSE_pctChg.pdf", width=6, height=5)


#--- Relative accuracy
relAcc.df <- comp %>% filter(!is.na(obs.lo)) %>%
  group_by(b, mtn, effort, era, spp) %>%
  summarise(obs.lo.relPrec=sd(obs.lo)/mean(obs.lo)*100,
            mod.lo.relPrec=sd(mod.lo)/mean(mod.lo)*100,
            obs.hi.relPrec=sd(obs.hi)/mean(obs.hi)*100,
            mod.hi.relPrec=sd(mod.hi)/mean(mod.hi)*100,
            obs.lo.relBias=mean((obs.lo-true.lo)/mean(true.lo))*100,
            mod.lo.relBias=mean((mod.lo-true.lo)/mean(true.lo))*100,
            obs.hi.relBias=mean((obs.hi-true.hi)/mean(true.hi))*100,
            mod.hi.relBias=mean((mod.hi-true.hi)/mean(true.hi))*100) %>%
  mutate(lo.diffPrec=mod.lo.relPrec - obs.lo.relPrec,
         hi.diffPrec=mod.hi.relPrec - obs.hi.relPrec,
         lo.diffBias=mod.lo.relBias - obs.lo.relBias,
         hi.diffBias=mod.hi.relBias - obs.hi.relBias) %>%
  mutate(obs.lo.relAcc=abs(obs.lo.relPrec) + obs.lo.relBias^2,
         mod.lo.relAcc=abs(mod.lo.relPrec) + mod.lo.relBias^2,
         obs.hi.relAcc=abs(obs.hi.relPrec) + obs.hi.relBias^2,
         mod.hi.relAcc=abs(mod.hi.relPrec) + mod.hi.relBias^2) %>%
  mutate(lo.relAccDiff=mod.lo.relAcc-obs.lo.relAcc,
         hi.relAccDiff=mod.hi.relAcc-obs.hi.relAcc) 
relAcc.df %>%
  select(b, mtn, effort, era, contains("relAccDiff")) %>%
  pivot_longer(5:6, names_to="Estimate", values_to="relAccDiff") %>%
  mutate(Boundary=str_sub(Estimate, 1, 2)) %>%
  mutate(Boundary=factor(Boundary, levels=c("lo", "hi"),
                         labels=c("Lower boundary", "Upper boundary")), 
         b=factor(b, levels=c("25", "50", "100", "200")),
         effort=factor(effort, levels=c(0.5, 1, 2), 
                       labels=c("50%", "100%", "200%")),
         era=factor(era, levels=c("H", "C"), 
                    labels=c("Historical", "Contemporary"))) %>%
  group_by(mtn, effort, era, Boundary, b) %>%
  summarise(mn=mean(relAccDiff, na.rm=T), 
            se=sd(relAccDiff, na.rm=T)/sqrt(n())) %>%
  ggplot(aes(x=effort, y=mn, ymin=mn-2*se, ymax=mn+2*se, colour=b, shape=mtn)) + 
  geom_hline(yintercept=0, linetype=2) + 
  geom_point(position=position_dodge(width=0.25), size=1.5) + 
  geom_linerange(position=position_dodge(width=0.25)) + 
  facet_grid(era~Boundary) +
  scale_shape_manual("Mountain\nRange", values=c(1,5)) +
  scale_colour_brewer("Bin size", type="qual", palette=2) +
  labs(x="Sampling effort relative to empirical", 
       y="Change in relative accuracy relative to empirical (mean ± 2 SE)") +
  theme_bw()
ggsave("figs/sims/RelAcc_pctChg.pdf", width=6, height=5)




