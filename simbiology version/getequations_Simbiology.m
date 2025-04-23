% Run script to print out all model equations and parameters
clear all;
% Load model
model = copyobj(sbioloadproject("CaRAS.sbproj").m1);
getequations(model) % print equations