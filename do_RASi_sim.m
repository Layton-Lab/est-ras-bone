function yvals = do_RASi_sim(do_ARB, do_ACEi, pct_ARB, pct_ACEi,...
                                age_ARB, age_ACEi, do_EST, tvals)
% outputs simulation results given settings
%    do_ARB: do ARB simulation (0 - false, 1 - true)
%    do_ACEi: do ACEi simulation (0 - false, 1 - true)
%    pct_ARB: percent inhibition by ARB 
%    pct_ACEi: percent inhibition by ACEi
%    age_ARB: age to start ARB
%    age_ACEi: age to start ACEi
%    do_EST: do estrogen decline (0 - false, 1 - true)
%    tvals: time values to integrate over
%% set initial conditions 
load("IC/2025-04-18_ICfinal3.mat","IC") % loads IC

%% set parameters
p = set_params;
[params, ~] = pars2vector(p,0);

%% Run simulations
% time span
t0 = 20*24*365; % Age 20 
tf = 80*24*365; % Age 80
tspan = [t0,tf];

% ode solver options
opts_ode = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9, 'MaxStep', 10);


[t, y] = ode15s(@(t,y) mod_eqns(t,y,params,...
                            'do_EST', do_EST,...
                            'do_ACEi', [do_ACEi,age_ACEi,pct_ACEi],...
                            'do_ARB', [do_ARB,age_ARB,pct_ARB]),...
                            tspan, IC, opts_ode);
tyrs = t/(365*24);

% Interpolate results over t-values
yvals = interp1(tyrs, y, tvals, 'linear');
end