# McCain et al. 2020
# Unusually large upward shifts in cold-adapted, montane mammals as temperature warms
# 
# Main Model


############--------------------------------------------------------------------
## Description
##-------------------------------------

# This script runs the hierarchical Bayesian undersampling model, partitioning
# the data into four subsets (2 mtn ranges x 2 eras). Note that the full model
# with adequate burn-in takes time to run








############--------------------------------------------------------------------
## Set up
##-------------------------------------

library(rjags); library(dclone); library(ggmcmc); library(tidyverse)

sets <- expand.grid(era=c("H", "C"),
                    mtn=c("FR", "SJ"))
sets.out <- vector("list", nrow(sets))








############--------------------------------------------------------------------
## Run model for each mtn x era combination
##-------------------------------------

for(i in 1:nrow(sets)) {
  
  #--- subset mountain, era
  era <- sets$era[i]
  mtn <- sets$mtn[i]  # mountain range
  
  
  #--- read in data
  # species information: filter by trapping method
  sp.i <- read_csv("data/Mamm_Summary_Data.csv") 
  pDet <- read_csv("data/prDet_processed.csv") %>%
    mutate(Abbrev=sp.i$Abbrev[match(Species, sp.i$Species)])
  # historic data: rename columns, select species
  Y.df <- read_csv(paste0("data/Mamm_", mtn, "_50m.csv")) %>%
    select("Elevation", contains(paste0("_", era))) %>%
    setNames(., str_remove(names(.), paste0("_", era))) %>%
    select("Elevation", one_of(sp.i$Abbrev))
  # interpolated ranges
  Y.rng <- map_dfc(2:ncol(Y.df), ~range(which(Y.df[,.]>0))) %>% as.matrix
  
  
  #--- reshape data for JAGS
  Y <- as.matrix(Y.df[,-1])  # remove Elevation column
  n.el <- nrow(Y)
  n.spp <- ncol(Y)
  Ytot <- rowSums(Y)
  spAbund <- colSums(Y)
  delta_shp <- as.matrix(pDet[match(names(spAbund), pDet$Abbrev),9:10])
  # calculate interpPatchy & binsAway for each species
  interpPatchy <- rep(0, n.spp)
  binsAway <- interpRng <- matrix(0, nrow=n.el, ncol=n.spp)
  for(s in 1:n.spp) {
    if(is.infinite(Y.rng[1,s])) {
      interpPatchy[s] <- 1
      binsAway[,s] <- n.el
    } else {
      interpRng.s <- Y.rng[1,s]:Y.rng[2,s]
      interpRng[,s] <- 1:n.el %in% interpRng.s
      interpPatchy[s] <- sum(Y[interpRng.s,s]==0)/length(interpRng.s)
      for(j in (1:n.el)[-interpRng.s]) {
        binsAway[j,s] <- min(abs(j - Y.rng[,s]))
      }
    }
  }
  
  
  #--- store data 
  jags_d <- list(J=n.spp, 
                 n.el=n.el, 
                 y=Y, 
                 Y=Ytot, 
                 delta_shp=delta_shp, 
                 delta_dat=delta_shp[,1]/rowSums(delta_shp),
                 LAMBDA=spAbund,
                 interpPatchy=c(scale(interpPatchy)),
                 distAway=matrix(scale(c(binsAway*50)), ncol=n.spp))
  obs.df <- data.frame(bin=as.character(rep(1:n.el, times=n.spp)),
                       Elevation=rep(Y.df$Elevation, times=n.spp),
                       spp=as.character(rep(1:n.spp, each=n.el)),
                       Abbrev=rep(names(spAbund), each=n.el),
                       Y=c(Y),
                       Ytot=rep(Ytot, times=n.spp),
                       binsAway=c(binsAway),
                       distAway=c(binsAway)*50,
                       distAway_sc=c(scale(c(binsAway)*50)),
                       interpPatchy=rep(interpPatchy, each=n.el),
                       spAbund=rep(spAbund, each=n.el),
                       delta=rep(delta_shp[,1]/rowSums(delta_shp), each=n.el),
                       era=era,
                       mtn=mtn,
                       set=paste0(mtn, "_", era))
  write_csv(obs.df, paste0("out/obs/obs_", i, ".csv"))
  
  
  #--- Run JAGS model
  pars <- c("lambda", "Z")
  cl <- makeCluster(3)
  out <- jags.parfit(cl=cl, data=jags_d, params=pars, 
                     model="code/models/1_hb_undersampling.txt",
                     inits=list(Z=matrix(1, n.el, n.spp)),
                     n.chains=3, n.adapt=50000, n.update=20000, 
                     n.iter=20000, n.thin=50)
  stopCluster(cl)
  
  
  #--- Summarize output
  sets.out[[i]] <- ggs(out, "Z") %>%
    mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
             str_remove("Z\\["),
           spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
             str_remove("\\]")) %>%
    group_by(bin, spp) %>%
    summarise(mn_prPres=mean(value)) %>%
    full_join(obs.df, c("bin", "spp")) 
  
}


#--- store output
write_csv(do.call('rbind', sets.out), "out/all_50m_bins.csv")




