# BN2MF for environmental epidemiology  

Bayesian non-parametric non-negative matrix factorization for environmental epidemiology.

## Roadmap

1. [sims](sims) has:
    * Script to create simulations.
    * Subfolder with csv files.
    * Subfolder with sims for bootstrap example.
    * `sim_sep.RDA` has all sims as a large nested dataframe.
    * `bs_sample.RDA` has 600 sim subsample for bootstrapping.
2. [main](main) has:
    * [bn2mf](bn2mf) subfolder with:
        * `run_bn2mf.m` to run on all sims.
        * `bn2mf_metrics.R` to calculate error metrics and coverage.
        * `combine_bn2mf.R` to combine metrics.
        * `output` subfolder with results.
    * [other_models](other_models) subfolder with:
        * `run_other_models` to run on all sims.
        * `other_model_metrics.R` to calculate error metrics.
        * `combine_other_models.R` to combine metrics.
        * `output` subfolder with results.
    * `make_figures.R` script combines all main results.
3. [bootstrap](bootstrap) has:
    * `bootstrap_bn2mf.m` to run on subsample.
    * `bootstrap_combo.R` to combine and reorder bootstrap results.
    * `bootstrap_coverage.R` to calculate coverage of true scores.
    * `ci_example.R` to vizualize results for one bootstrapped example.
    * `bootstrap_nmf.R` to bootstrap Poisson NMF.
    * `output` subfolder with results.
4. [functions](functions) has functions I've written, see code for comments.
5. [figures](figures) has figures for manuscript.
6. [misc](misc) has miscellaneous scripts.
    * Related, but not necessary to run analysis.
    * Feel free to ignore.
