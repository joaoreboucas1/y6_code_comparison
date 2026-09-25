# Code comparison notes
## Sep 25, 2026:
- Currently disabling the ell_prefactor2 = l*(l-1)*(l+1)*(l+2)/(l+0.5)^4 in Cocoa (see eqs. 78 and 79 in https://arxiv.org/pdf/1812.05995)
- Disabling the ell_prefactor2 gets rid of large-angle drifts between Cocoa and Cosmolike/Cosmosis, this is what i'm adopting from now on
- Cosmosis uses `alpha` and Cosmolike uses `B_MAG = 2*(alpha - 1)`, this is fixed now
- Lighthouse has a bug where the shear calibration parameters are set but not taken into account in the calculation
- Cosmolike and Cocoa agree on cosmic shear, and dchi2 between cosmosis and cosmolike is ~0.36. now we want to invetigate magnification bias, so we are turning it off from both analysis. This means alpha = 1 and BMAG = 0.
- Removing BMAG, dchi2 in clustering between Cocoa and Cosmolike went down from 1.1 to 0.23!
- Removing a bmag_ell_prefactor in Cocoa
- Changing line 2718 in cosmo2D.c in Cocoa to integrate until a = 1 instead of amax