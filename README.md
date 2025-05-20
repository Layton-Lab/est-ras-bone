# RAS, estrogen, calcium, and bone
This code is for the manuscript [Stadt and Layton, "A physiology-based mathematical model of the renin-angiotensin system, bone remodeling, and calcium homeostasis: Effects of estrogen and renin-angiotensin system inhibitors", submitted](https://www.biorxiv.org/content/10.1101/2025.05.01.651663v2.abstract)
 The manuscript is under review.


## Key files
**driver_ANGandPTHinf.m** run this script to conduct the PTH and Ang II infusion experiments

**driver_multi_RASi.m** run this script to run simulations of age-related estrogen decline with and without RAS inhibitors of varying efficacy

**driver_vary_age.m** run this script to run simulations of age-related estrogen decline with RASi starting at varied ages

**driverSS.m** run this script to compute the model steady state

**./simbiology version** this directory contains a Simbiology version of the model and scripts for running Sobol sensitivity analysis
