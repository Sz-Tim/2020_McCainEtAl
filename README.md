### McCain CM, SRB King, TM Szewczyk. 2020. Unusually large upward shifts in cold-adapted, montane mammals as temperature warms. *Ecology*.

---  

Code repository accompanying the manuscript.


### Table of Contents  
- **data/** contains data files, including raw files and files produced by scripts in this repository.
  - Mamm_Summary_Data.csv: Elevational minimum, maximum, and number of detections for each species in each moutain range and era  
  - Mamm_FR_50m.csv: Detection data for the Front Range at 50m elevational bins  
  - Mamm_SJ_50m.csv: Detection data for the San Juans at 50m elevational bins  
  - prDet_processed.csv: Prior distributions for species' detection probabilities  
  - sample_els_C.csv: County, elevation, and species for contemporary detections  
  - sample_els_H.csv: County, elevation, and species for historic detections  
  - occupancy_mtn.csv: Summaries from repeatedly sampled sites for comparison with occupancy models  
  - CO_SiteInfo.txt: Basic information for repeatedly sampled 2010-2012 sites (Appendix S2)  
  - orig/: .xlsx files with a sheet for each species
  - sites/: .txt files with data for repeatedly sampled sites  
- **code/** contains R scripts for processing data and running models  
  - 00_fn.R: Functions for simulations  
  - 1_dataProcessing.R: Code for munging raw data  
  - 2_mainModel.R: Code for running the undersampling model for each era x mtn  
  - 3_priorPredictions.R: Code for prior predictive checks  
  - 4_sensitivityAnalysis.R: Code for simulating communities & assessing impact of elevational bin size and sampling effort  
  - 5_occupancyComparison.R: Code for running the undersampling model and an occupancy model for comparison with only the repeatedly sampled sites  
  - models/: JAGS models
- **out/** output produced by running the scripts in code/  
- **figs/** select figures produced by running the scripts in code/  



