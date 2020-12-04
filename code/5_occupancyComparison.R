# McCain et al. 2020
# Unusually large upward shifts in cold-adapted, montane mammals as temperature warms
# 
# Occupancy Model Comparison


############--------------------------------------------------------------------
## Description
##-------------------------------------

# This script runs the hierarchical Bayesian undersampling model and an 
# elevational occupancy model for comparison, using the subset of the 
# contemporary data that derives from discrete sites that were sampled across
# five nights (Appendix S1: Contemporary Mammal Sampling)








############--------------------------------------------------------------------
## Set up
##-------------------------------------

library(tidyverse); library(rjags); library(ggmcmc); 
library(dclone); library(HDInterval)

sp.i <- read_csv("data/Mamm_Summary_Data.csv")
pDet <- read_csv("data/prDet_processed.csv") %>%
  mutate(Abbrev=sp.i$Abbrev[match(Species, sp.i$Species)])
occ.df <- read_csv("data/occupancy_mtn.csv")








############--------------------------------------------------------------------
## Run models for each mtn
##-------------------------------------

gg.occ <- gg.hbu <- vector("list", 2)

for(i in 1:2) {
  
  mtn_i <- unique(occ.df$mtn)[i]
  occ_i <- filter(occ.df, mtn==mtn_i)
  
  # all elevations on gradient binned to 50m
  el_i <- seq(min(occ_i$ElBin), max(occ_i$ElBin), by=50)
  n.el <- length(el_i)
  J_i <- n_distinct(occ_i$Species)
  y_i <- occ_i %>% select(ElBin, Species, nDet) %>%
    pivot_wider(names_from=Species, values_from=nDet, values_fill=0) %>%
    bind_rows(tibble(ElBin=el_i[!el_i %in% unique(occ_i$ElBin)])) %>%
    replace(is.na(.), 0) %>% arrange(ElBin)
  k_i <- occ_i %>% group_by(ElBin) %>% 
    summarise(k=max(k)) %>% ungroup %>% 
    bind_rows(tibble(ElBin=el_i[!el_i %in% unique(occ_i$ElBin)])) %>%
    replace(is.na(.), 0) %>% arrange(ElBin)
  n_i <- occ_i %>% select(ElBin, Species, nInd) %>%
    pivot_wider(names_from=Species, values_from=nInd, values_fill=0) %>% 
    bind_rows(tibble(ElBin=el_i[!el_i %in% unique(occ_i$ElBin)])) %>%
    replace(is.na(.), 0) %>% arrange(ElBin)
  
  # detection priors
  delta_shp <- as.matrix(pDet[match(colnames(y_i)[-1], pDet$Species),9:10])
  
  # patchiness & distance from boundaries
  interpPatchy <- rep(0, J_i)
  binsAway <- interpRng <- distAway <- matrix(0, nrow=n.el, ncol=J_i)
  Y.rng <- apply(y_i[,-1], 2, function(x) range(which(x>0)))
  for(s in 1:J_i) {
    interpRng.s <- Y.rng[1,s]:Y.rng[2,s]
    interpRng[,s] <- 1:n.el %in% interpRng.s
    interpPatchy[s] <- sum(y_i[interpRng.s,s+1]==0)/length(interpRng.s)
    for(j in (1:n.el)[-interpRng.s]) {
      binsAway[j,s] <- min(abs(j - Y.rng[,s]))
      distAway[j,s] <- min(abs(el_i[j] - el_i[Y.rng[,s]]))
    }
  }
  
  
  
  
  jags_d <- list(J=J_i, 
                 n.el=length(el_i), 
                 el=c(scale(el_i)),
                 el2=c(scale(el_i))^2,
                 k=k_i$k,
                 y_det=as.matrix(y_i[,-1]), 
                 y=as.matrix(n_i[,-1]),
                 Y=rowSums(as.matrix(n_i[,-1])), 
                 delta_shp=delta_shp,
                 LAMBDA=colSums(as.matrix(n_i[,-1])),
                 interpPatchy=c(scale(interpPatchy)),
                 distAway=matrix(scale(c(distAway)), ncol=J_i))
  
  obs.df <- data.frame(bin=rep(1:jags_d$n.el, times=jags_d$J),
                       Elevation=rep(el_i, times=jags_d$J),
                       spp=rep(1:jags_d$J, each=jags_d$n.el),
                       binsAway=c(binsAway),
                       Y=c(jags_d$y),
                       Ytot=rep(jags_d$Y, times=jags_d$J),
                       mtn=mtn_i, 
                       Species=rep(colnames(y_i)[-1], each=jags_d$n.el))
  
  
  
  
  #--- occupancy model
  cl <- makeCluster(3)
  out.occ <- jags.parfit(cl=cl, data=jags_d, params=c("Z", "psi"),
                         model="code/models/2_occupancy.txt", 
                         inits=list(Z=matrix(1, jags_d$n.el, jags_d$J)),
                         n.chains=3, n.adapt=50000, n.update=20000, 
                         n.iter=20000, n.thin=50)
  stopCluster(cl)
  
  gg.occ[[i]] <- full_join(
    ggs(out.occ, "psi") %>%
      mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
               str_remove("psi\\[") %>% as.numeric,
             spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
               str_remove("\\]") %>% as.numeric) %>%
      group_by(bin, spp) %>% 
      summarise(psi_mn=mean(value)), 
    ggs(out.occ, "Z") %>%
      mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
               str_remove("Z\\[") %>% as.numeric,
             spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
               str_remove("\\]") %>% as.numeric) %>%
      group_by(bin, spp) %>%
      summarise(Z_mn=mean(value)), 
    by=c("bin", "spp")) %>%
    full_join(., obs.df, c("bin", "spp")) 
  
  
  
  
  #--- hierarchical bayesian undersampling model
  cl <- makeCluster(3)
  out.hbu <- jags.parfit(cl=cl, data=jags_d, params=c("Z", "psi"), 
                         model="code/models/1_hb_undersampling.txt",
                         inits=list(Z=matrix(1, jags_d$n.el, jags_d$J)),
                         n.chains=3, n.adapt=500, n.update=200, 
                         n.iter=200, n.thin=5)
  stopCluster(cl)
 
  gg.hbu[[i]] <- full_join(
    ggs(out.hbu, "psi") %>%
      mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
               str_remove("psi\\[") %>% as.numeric,
             spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
               str_remove("\\]") %>% as.numeric) %>%
      group_by(bin, spp) %>% 
      summarise(psi_mn=mean(value)),
    ggs(out.hbu, "Z") %>%
      mutate(bin=str_split_fixed(Parameter, ",", 2)[,1] %>%
               str_remove("Z\\[") %>% as.numeric,
             spp=str_split_fixed(Parameter, ",", 2)[,2] %>%
               str_remove("\\]") %>% as.numeric) %>%
      group_by(bin, spp) %>%
      summarise(Z_mn=mean(value)), 
    by=c("bin", "spp")) %>%
    full_join(., obs.df, c("bin", "spp")) 
  
  
}








############--------------------------------------------------------------------
## Summarise and plot
##-------------------------------------

gg.comp <- list("FR"=rbind(mutate(gg.occ[[1]], mod="occ"), 
                           mutate(gg.hbu[[1]], mod="hbu")),
                "SJ"=rbind(mutate(gg.occ[[2]], mod="occ"), 
                           mutate(gg.hbu[[2]], mod="hbu")))

gg.comp.df <- do.call("rbind", gg.comp) %>% 
  mutate(Detection=ifelse(Y>0, "Detection", "Non-detection"),
         Species_=str_replace(Species, " ", "\n"), 
         Species_=factor(Species_, levels=rev(sort(unique(Species_)))))

p.psi <- ggplot(gg.comp.df, aes(Elevation, y=Species_)) +
  geom_vline(data=filter(gg.comp.df, Ytot>0), aes(xintercept=Elevation), 
             colour="gray60", size=0.25, linetype=3) +
  ggridges::geom_ridgeline(aes(height=psi_mn, fill=mod, colour=mod), 
                           alpha=0.5, scale=0.5) +
  geom_point(data=filter(gg.comp.df, Detection=="Detection"), 
             aes(shape=Detection), 
             colour="black", size=0.75, position=position_nudge(y=0.5)) +
  scale_colour_manual("Model", values=c("red", "cadetblue"), 
                      labels=c("Hierarchical\nUndersampling",
                               "Occupancy")) +
  scale_fill_manual("Model", values=c("red", "cadetblue"), 
                    labels=c("Hierarchical\nUndersampling",
                             "Occupancy")) +
  scale_shape_manual("", values=c(19,4)) +
  facet_wrap(~mtn, scales="free_x") +
  theme_bw() +
  theme(panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(), 
        panel.grid.major.y=element_line(colour="gray30"),
        axis.text.y=element_text(vjust=0, lineheight=0.7, size=11),
        axis.text.x=element_text(size=11), 
        axis.title.x=element_text(size=12),
        plot.title=element_text(size=14), 
        strip.text=element_text(size=12),
        legend.text=element_text(size=11), 
        legend.title=element_text(size=12)) + 
  labs(x="Elevation (m)", y="", 
       title=expression(A.~Probability~of~presence:~psi))

p.Z <- ggplot(gg.comp.df, aes(Elevation, y=Species_)) +
  geom_vline(data=filter(gg.comp.df, Ytot>0), aes(xintercept=Elevation), 
             colour="gray60", size=0.25, linetype=3) +
  ggridges::geom_ridgeline(aes(height=Z_mn, fill=mod, colour=mod), 
                           alpha=0.5, scale=0.5) +
  geom_point(data=filter(gg.comp.df, Detection=="Detection"), 
             aes(shape=Detection), 
             colour="black", size=0.75, position=position_nudge(y=0.5)) +
  scale_colour_manual("Model", values=c("red", "cadetblue"), 
                      labels=c("Hierarchical\nUndersampling",
                               "Occupancy")) +
  scale_fill_manual("Model", values=c("red", "cadetblue"), 
                    labels=c("Hierarchical\nUndersampling",
                             "Occupancy")) +
  scale_shape_manual("", values=c(19,4)) +
  facet_wrap(~mtn, scales="free_x") +
  theme_bw() +
  theme(panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(), 
        panel.grid.major.y=element_line(colour="gray30"),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=11), 
        axis.title.x=element_text(size=12),
        plot.title=element_text(size=14), 
        strip.text=element_text(size=12),
        legend.text=element_text(size=11), 
        legend.title=element_text(size=12)) + 
  labs(x="Elevation (m)", y="", 
       title=expression(B.~Latent~predicted~presence:~Z))

ggpubr::ggarrange(p.psi, p.Z, widths=c(1, 0.8),
                  common.legend=T, legend="bottom") %>%
  ggsave("figs/comp_Occupancy.jpg", ., height=12, width=10, dpi=300, units="in")

