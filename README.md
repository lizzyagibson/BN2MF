# BN2MF for environmental epidemiology  

Bayesian non-parametric non-negative matrix factorization for environmental epidemiology.

## roadmap

*ignore `.m` files.*  
    
1. [sims](sims) has:  
    * `sim_grid.R` script to create simulations. (1st)  
    * `sim_sep.RDA` with all sims as a large nested dataframe.  
        * `sim_ids.RDA` with just id variables.  
    * `bs_sample.RDA` with 600 sim subsample for bootstrapping.  
            * `bs_ids.csv` with just id variables.  
    * `csvs` subfolder with one example simulated dataset.  
2. [functions](functions) has functions I've written, see code for comments.  
    * `compare_functions.R` has functions that I use throughout analysis. (2nd)  
    * **ignore `fig_set.R`.**  
3. [main](main) has:  
    * [bn2mf](bn2mf) subfolder with:  
        * `run_bn2mf.m` to run on all sims.  
        * `bn2mf_metrics.R` to calculate error metrics and coverage. (3rd)  
        * `output` subfolder with results.  
    * [other_models](other_models) subfolder with:  
        * `run_other_models` to run on all sims. (4th)  
        * `other_model_metrics.R` to calculate error metrics. (5th)  
        * `output` subfolder with results.  
    * `combine_main.R` to combine bn2mf and other model metrics. (6th)  
    * `make_figures.R` script combines all main results. (7th)  
4. [bootstrap](bootstrap) has:  
    * three files for BN2MF bootstrap:  
        * `bootstrap_bn2mf.m` to run on subsample.  
        * `bootstrap_combo.R` to combine and reorder bootstrap results. (8th)  
        * `bootstrap_coverage.R` to calculate coverage of true scores. (9th)  
    * three files for NMF bootstrap (**ignore**):  
        * `bootstrap_nmf.R` to bootstrap Poisson NMF.  
        * `nmf_bootstrap_combo.R` to combine and reorder bootstrap results.  
        * `nmf_bootstrap_coverage.R` to calculate coverage of true scores.  
    * `vci_bs_coverage.R` saves VCI metrics to compare with bootstraps. (10th)  
    * `compare_bootstrap.R` summarizes both bootstraps and VCI. (11th)  
    * For 8-10, you can't run these, but please check for logic.  
    * `output` subfolder with results.  
5. [misc](misc) has miscellaneous scripts (**ignore**). 
    * Related, but not necessary to run analysis.

