# RAS, estrogen, calcium, and bone
This code is for the manuscript ``A physiology-based mathematical model of the renin-angiotensin system, bone remodeling, and calcium homeostasis: Effects of estrogen and renin-angiotensin system inhibitors''
by Melissa Stadt and Anita Layton. The manuscript is under review.


## Key files
**getequations_Simbiology** prints out all model equations to command line

**driver_estdecline.m** runs model simulations with age-related estrogen decline

**setSobol_newpars.mlx** use this to set up a .mat file to use with "runSobol.m". save results under "SobolSet" and change file name in "runSobol.m" to run Sobol analysis

**runSobol.m** run Sobol analysis given file from "SobolSet"

**makeplot_sobolresults.m** make the plots from sobol results (NOTE: need to run "runSobol.m" first to get results)

