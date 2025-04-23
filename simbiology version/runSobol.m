% Run Sobol GSA
clearvars;

% Load Sobol GSA settings (set with setSobol.mlx)
load("./SobolSet/21-Apr-2025_SobolSetings_notes-newpars_final.mat")

% Run Sobol
fprintf("start Sobol \n")
tic;
sobolResults = sbiosobol(model, parnames, obs,...
                "Bounds", bounds,...
                "UseParallel", true,...
                "Accelerate", true,...
                "ShowWaitBar", true,...
                "OutputTimes", output_times);
elapsedTime = toc;
fprintf("Elapsed time is %.2f mins. \n", elapsedTime/60)

% Save results
fname = strcat(date, "_runSobol_CaRAS",...
    "_notes-", notes,...
    ".mat");
save(fname)

fprintf("done! \n results saved to \n %s \n", fname);