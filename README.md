# DES-Y6 Cosmolike-Cosmosis Code Comparison

This repository contains files to compare the Cosmolike, Cosmosis and Cocoa codes regarding the analysis of DES-Y6 data.

## Contents

- `cosmosis` contains the Cosmosis setup (`ini/`, `modules_desy6/`, `inc_ini_nzmode/`), the Y6 data file (`EXTENSION_DVv8.0_fixmnu.fits`), `parse_cov.py`, and the pipeline outputs for the LCDM and w0waCDM fiducial cosmologies (`lcdm_datavector_run/`, `w0wacdm_datavector_run/`). See `cosmosis/README.md` for the fiducial parameters.
- `cosmolike_lighthouse` contains model vectors generated with Cosmolike Lighthouse (`LCDM_lighthouse*.modelvector`)
- `cocoa` contains model vectors generated with Cocoa (`LCDM*.modelvector`), the evaluation yaml (`EXAMPLE_EVALUATE1.yaml`), and its outputs
- `comparing` contains the Cosmosis model vector in Cosmolike format (`COSMOSIS.modelvector`), the Y6 covariance and mask (`DESY6.cov`, `DESY6.mask`), and the comparison tools: `calc_chi2_and_plot.py` computes Delta chi2 (3x2pt, 2x2pt, ggl, wtheta, shear) and difference plots for each pair of codes, saved to `plots/`, and `plot.ipynb`/`utils.py` make the per-probe comparison plots

Notes on the comparison process are given in `NOTES.md`.

Cocoa users can also take a look at [cocoa_des_y6](https://github.com/joaoreboucas1/cocoa_des_y6) repository which implements DES-Y6 but is still a work in progress.