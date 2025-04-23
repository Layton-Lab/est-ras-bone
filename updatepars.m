% run this script to update the parameter list to put in "mod_eqns"
% writes parameter list to "newparamnames.txt"
clear all;

%% Set parameters
p = set_params();
% Change values here

[params, parnames] = pars2vector(p,1);