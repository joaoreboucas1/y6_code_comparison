# DES-Y6 Cosmolike-Cosmosis Code Comparison

This repository contains files to compare the Cosmolike, Cosmosis and Cocoa codes regarding the analysis of DES-Y6 data. The main results of the comparison are in `comparing/plots/`. A Claude write-up of the progress and results is in https://claude.ai/artifact/6YKN3GSrjzkfbYD9deDmyg.

For the comparison, we used:
- Modified version of Cocoa v4.11.2 where the variables `ell_bmag_prefactor` and `ell_prefactor2` are set to one; using Cocoa project [cocoa_des_y6](https://github.com/joaoreboucas1/cocoa_des_y6) which implements DES-Y6 but is still a work in progress;
- Cosmolike_core branch cluster_chto commit 1b861b1fbe73a2d62dfce3acf29970cff20f5e30 with a slight modification: in `init.c` and `init_basic.c`, we set `nuisance.oneplusz0_ia=1.62` and `nuisance.c1rhocrit_ia=0.013873073650776856`; [Lighthouse](git@github.com:CosmoLike/lighthouse.git) branch `des_y6_comparison` (mostly a cleanup of the `Ystatistics` branch);
- Cosmosis: DES version of cosmosis-standard-library commit 54b15ffb9dd1a71cf6a936c752bb9c8df143beeb.

## Contents

- `cosmosis` contains the Cosmosis setup (`ini/`, `modules_desy6/`, `inc_ini_nzmode/`), the Y6 data file (`EXTENSION_DVv8.0_fixmnu.fits`), `parse_cov.py`, and the pipeline outputs for the LCDM and w0waCDM fiducial cosmologies (`lcdm_datavector_run/`, `w0wacdm_datavector_run/`). See `cosmosis/README.md` for the fiducial parameters.
- `cosmolike_lighthouse` contains model vectors generated with Cosmolike Lighthouse (`LCDM_lighthouse*.modelvector`)
- `cocoa` contains model vectors generated with Cocoa (`LCDM*.modelvector`), the evaluation yaml (`EXAMPLE_EVALUATE1.yaml`), and its outputs
- `comparing` contains the Cosmosis model vector in Cosmolike format (`COSMOSIS.modelvector`), the Y6 covariance and mask (`DESY6.cov`, `DESY6.mask`), and the comparison tools: `calc_chi2_and_plot.py` computes Delta chi2 (3x2pt, 2x2pt, ggl, wtheta, shear) and difference plots for each pair of codes, saved to `plots/`, and `plot.ipynb`/`utils.py` make the per-probe comparison plots

Notes on the comparison process are given in `NOTES.md`.