# S2S_PostProcess

0_control_PLSRmdl.R 
- input for all other scripts
- sets domain and 'cut' number
- sets all directories used in scripts
- sets variables names
- cross validation parameters

1_PLSR_data_prep... 
- file which prepares the predictors for processing with PLSR

2_... 
- all files perform PLSR on hucs in some for; many of these were used for initial testing
- 'allHRU_PLSR_multiVar_allLeadsMons' run all huc/leads/mons in parallel (one file doesnt do CV)

3_...
- processes PLSR output
- some do basic processing, others do full analysis

4_... 
- plots to visualize results
