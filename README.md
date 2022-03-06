# plectranthus_barbatus
This code reproduces the results from "Mechanistic-statistical modeling reveals the pivotal role of human-mediated long-distance dispersal in plant invasion dynamics"

Guidelines to reproduce:
1) Download data files (**data_for_model.Rdata, spatial_cells.Rdata, validation_kit.Rdata**) our Zenodo repository: https://doi.org/10.5281/zenodo.6331676
2) Copy R scripts of this repository and data files in a common folder.
3) Optional: Run script **fit_model.R** which will load the data file **data_for_model.R**. If you run all chains sequentially as proposed by default in the script, it should take at least a week to be completed. You should obtain exactly the files **samples_pt21_seed32.Rdata** etc. which are also provided in this repository for a quick reproduction.  
4) Run script **output_Figures.R** which will load all data files, create and save the Figures of our manuscript as .png files.
