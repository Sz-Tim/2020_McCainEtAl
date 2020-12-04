# McCain et al. 2020
# Unusually large upward shifts in cold-adapted, montane mammals as temperature warms
# 
# Data Processing


############--------------------------------------------------------------------
## Description
##-------------------------------------

# This script processes raw data for use in the models.
# 1) Aggregate datasets
#      Raw detection data is aggregated for historic and contemporaneous 
#      datasets and saved to 'data/sample_els_H.csv' and 'data/sample_els_C.csv'
#      where each contains a row for every observation and columns for the 
#      county (in Colorado, USA), elevation (m), and species abbreviation
# 2)  Estimate detection distributions
#      Detection probabilities from each repeatedly sampled site in the 
#      contemporaneous dataset are loaded, and beta distributions are estimated
#      for each species, which then serve as the prior distributions in the
#      hierarchical undersampling model. We use fitdistrplus::fitdist(), which 
#      uses same parameterization as stats::Beta (shape1, shape2). The dataframe
#      is saved to 'data/prDet_processed.csv' and contains a row for each 
#      species and columns describing the detection probability distribution
# 3) Occupancy model data subset
#      Loads data from the 32 repeatedly sampled sites in the contemporaneous
#      dataset and formats them for the occupancy model comparison. The 
#      dataframe is saved to 'data/occupancy_mtn.csv' and contains a row for 
#      each species detected in each 50m elevational bin in each mountain
#      range, and columns for the number of detections, number of trap-nights, 
#      and number of individuals detected








library(fitdistrplus); library(tidyverse); library(readxl)


############--------------------------------------------------------------------
## Aggregate datasets
##-------------------------------------

f.path <- dir("data/orig", "*ata.xlsx", full.names=T)
H.orig.ls <- C.orig.ls <- vector("list", length(f.path))

for(f in f.path) {
  f.sheets <- excel_sheets(f) %>%
    grep(pattern="Sampling", value=T, invert=T) %>% 
    grep(pattern="UNK", value=T, invert=T)
  H.orig.ls[[f]] <- map_dfr(f.sheets, 
                         ~read_xlsx(f, .x, range=cell_cols("A:B"),
                                    col_types=c("text", "numeric")) %>%
                           setNames(c("county", "el")) %>% mutate(sp=.x))
  C.orig.ls[[f]] <- map_dfr(f.sheets, 
                         ~read_xlsx(f, .x, range=cell_cols("I:J"),
                                    col_types=c("text", "numeric")) %>%
                           setNames(c("county", "el")) %>% mutate(sp=.x))
}

H.orig.df <- do.call("rbind", H.orig.ls)
C.orig.df <- do.call("rbind", C.orig.ls)

write_csv(H.orig.df, "data/sample_els_H.csv")
write_csv(C.orig.df, "data/sample_els_C.csv")








############--------------------------------------------------------------------
## Estimate detection distributions
##-------------------------------------

# probability of detection
# raw data: rows = species (44), columns = sites (32)
prDet.raw <- read_xlsx("data/orig/Prob_Detection.xlsx", 1) %>%
  mutate(Genus=str_split_fixed(Species, " ", 2)[,1]) 

# long-form, with only species that were detected (24) at the 32 sites
prDet.df <- prDet.raw %>%
  select(Order, Family, Genus, Species, 5:36) %>% 
  pivot_longer(5:36, names_to="Site", values_to="prDet", values_drop_na=T) 

# species-level: for species with ≥2 values
prDet.sum <- prDet.df %>%  group_by(Genus, Species) %>%
  summarise(tryGenus=n()<2, useFam=F, mn=mean(prDet), sd=sd(prDet), 
            shp1=ifelse(tryGenus, NA,
                        fitdist(prDet, "beta",
                                start=list(shape1=3, shape2=5))$estimate[1]),
            shp2=ifelse(tryGenus, NA,
                        fitdist(prDet, "beta",
                                start=list(shape1=3, shape2=5))$estimate[2]))

# genus-level: if <2 values wihin species, but ≥2 values within genus
prDet.genus <- prDet.df %>%  group_by(Family, Genus) %>%
  summarise(useFam=n()<2, mn=mean(prDet), sd=sd(prDet),
            shp1=ifelse(useFam, NA,
                        fitdist(prDet, "beta",
                                start=list(shape1=3, shape2=5))$estimate[1]),
            shp2=ifelse(useFam, NA,
                        fitdist(prDet, "beta",
                                start=list(shape1=3, shape2=5))$estimate[2]))

# family-level: if <2 values wihin genus
prDet.fam <- prDet.df %>%  group_by(Family) %>%
  summarise(useOrder=n()<2, mn=mean(prDet), sd=sd(prDet),
            shp1=ifelse(useOrder, NA,
                        fitdist(prDet, "beta",
                                start=list(shape1=3, shape2=5))$estimate[1]),
            shp2=ifelse(useOrder, NA,
                        fitdist(prDet, "beta",
                                start=list(shape1=3, shape2=5))$estimate[2]))

# add all taxa (i.e., 20 undetected species), info
prDet.sp <- prDet.raw %>% select(Order, Family, Genus, Species) %>%
  left_join(., filter(prDet.sum, !is.na(shp1)), by=c("Species", "Genus")) %>%
  mutate(tryGenus=replace(tryGenus, is.na(tryGenus), TRUE))
prDet.genus <- prDet.raw %>% group_by(Genus, Family, Order) %>% summarise() %>%
  left_join(., prDet.genus %>% filter(!is.na(shp1)), 
            by=c("Family", "Genus")) %>%
  mutate(useFam=replace(useFam, is.na(useFam), TRUE))

# merge dataframes: rows = species (44) with Beta shape parameters
prDet.full <- rbind(
  prDet.sp %>% filter(!tryGenus), # species values
  left_join(prDet.sp %>% 
              filter(tryGenus) %>% 
              select(-useFam, -mn, -sd, -shp1, -shp2), 
            rbind(prDet.genus %>%  # genus values 
                    filter(!useFam) %>% ungroup,
                  left_join(prDet.genus %>% 
                              filter(useFam) %>% ungroup %>%
                              select(-mn, -sd, -shp1, -shp2),
                            prDet.fam %>% # family values
                              select(-useOrder), 
                            by="Family")
            ),
            by=c("Order", "Family", "Genus")
  )
)

with(prDet.full, sum(!tryGenus)) # 18 spp with species-level distributions
with(prDet.full, sum(tryGenus & !useFam)) # 15 spp with genus-level distributions
with(prDet.full, sum(useFam)) # 11 spp with family-level distributions

write_csv(prDet.full, "data/prDet_processed.csv")

# visualize Beta distributions
cb.palette <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
for(f in 1:n_distinct(prDet.full$Family)) {
  fam <- unique(prDet.full$Family)[f]
  gen <- unique(filter(prDet.full, Family==fam & tryGenus)$Genus)
  gen <- setNames(seq_along(gen), gen)
  gray.palette <- paste0("gray", round(seq(30, 80, length.out=length(gen))))
  pDet.f <- filter(prDet.full, Family==fam) %>%
    mutate(lty=rep(c(1,5), each=length(cb.palette))[1:n()], 
           lty=ifelse(tryGenus, 2, lty),
           lty=ifelse(useFam, 3, lty),
           col=rep(cb.palette, 10)[1:n()], 
           col=ifelse(tryGenus, gray.palette[gen[Genus]], col),
           col=ifelse(useFam, "black", col),
           lwd=ifelse(tryGenus, 2, 1.5),
           lwd=ifelse(useFam, 1.5, lwd)) %>%
    arrange(Species)
  pdf(paste0("figs/pDet/pDet_priors_", fam, ".pdf"), width=7, height=7)
  plot(NA, NA, xlim=c(0,1), ylim=c(0,20), 
       xlab="Detection probability", ylab="Density", main=fam)
  walk(1:nrow(pDet.f), ~curve(dbeta(x, pDet.f$shp1[.], pDet.f$shp2[.]), 
                              from=0, to=1, add=T, col=pDet.f$col[.], 
                              lty=pDet.f$lty[.], lwd=pDet.f$lwd[.]))
  legend("topright", pDet.f$Species, col=pDet.f$col, lty=pDet.f$lty,
         bty="n", lwd=pDet.f$lwd)
  dev.off()
}








############--------------------------------------------------------------------
## Occupancy model data subset
##-------------------------------------

# original trap data
sp.i <- read_csv("data/Mamm_Summary_Data.csv") 
site_i <- read_csv("data/CO_SiteInfo.txt") 
trap.f <- dir("data/sites", ".txt")  # file for each site

# read in all files
trap.df <- map_dfr(setNames(trap.f, str_remove(str_sub(trap.f, 1, -5), " ")), 
                   ~read_csv(paste0("data/sites/", .x)) %>% 
                     select(Date, Species, `Tag/ Clip`) %>%
                     mutate(Tag=as.character(`Tag/ Clip`)) %>%
                     select(-`Tag/ Clip`), .id="Site") %>%
  # standardize species names across all files
  mutate(Date=lubridate::mdy(str_split_fixed(Date, " ", 2)[,1]),
         Species=str_remove(Species, "\\?"),
         Species=str_remove(Species, "\\(\\)"), 
         Species=str_remove(Species, "\\."),
         Species=str_trim(Species)) %>%
  mutate(Species_=case_when(grepl("gapperi", Species) ~ "Myodes gapperi",
                            Species=="Neotoma mexicanus" ~ "Neotoma mexicana",
                            Species=="Spermophilus lateralis" ~ "Callospermophilus lateralis",
                            Species=="Spermophilus variegatus" ~ "Otospermophilus variegatus",
                            Species=="Tamiascuius hudsonicus" ~ "Tamiasciurus hudsonicus")) %>%
  mutate(Species_=if_else(is.na(Species_), Species, Species_)) %>%
  select(-Species) %>% rename(Species=Species_) %>%
  filter(Species %in% sp.i$Species) %>%
  left_join(., site_i, by="Site") %>%
  mutate(ElBin=round(ElAvg/50)*50,
         mtn=case_when(Transect %in% c("Big Thompson", "Boulder") ~ "FR",
                       Transect %in% c("Dolores", "Lizardhead") ~ "SJ"))

# summarise number of detections, trap-nights, and individuals
trap_det.df <- trap.df %>% 
  group_by(mtn, Site, ElBin, Species, Date) %>% 
  summarise(nDet=n()) %>% 
  group_by(mtn, Site, ElBin, Species) %>%
  summarise(nDet=sum(nDet>0), k=5) %>% 
  group_by(mtn, ElBin, Species) %>%
  summarise(nDet=sum(nDet), k=sum(k)) %>%
  left_join(., 
            trap.df %>%
              group_by(mtn, ElBin, Species) %>%
              summarise(nInd=n_distinct(Tag)) %>% ungroup,
            by=c("mtn", "ElBin", "Species"))

write_csv(trap_det.df, "data/occupancy_mtn.csv")
