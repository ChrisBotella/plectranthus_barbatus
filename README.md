# plectranthus_barbatus
This code reproduces the results and Figures from "Dynamic species distribution modeling reveals the pivotal role of human-mediated long-distance dispersal in plant invasion"

Guidelines to reproduce:
1) Download data files (**data_for_model, data_full, validation_kit, toAnalyse_Validation.Rdata** and pre-fitted posterior samples if needed) from our Zenodo repository (V2): https://doi.org/10.5281/zenodo.6921965
2) Download R scripts from this repository in the same folder as the data files.
3) Optional: Run script **fit_model.R** which will load the data files **data_for_model.R** and **data_full** to fit the model on the full data. If you run all chains sequentially as proposed by default in the script, it should take at least a week to be completed. You should obtain exactly the Rdata files **samples_pt33_seed1**, **samples_pt33_seed1_logLik**, etc. which are also provided in the Zenodo repository for a quicker reproduction.  
4) Run script **convergence_check.R** which will load all posterior sample files, thin them (and save as "toAnalyse.Rdata"), produce trace plots and compute convergence criterions (see Appendix S3).
5) Run script **output_Figures.R** which will load all data files, create and save all the Figures of our manuscript as .png files.
