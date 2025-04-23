% Run simulations of estrogen decline with and without RASi
clearvars; % clear

% Load model
model = copyobj(sbioloadproject("CaRAS.sbproj").m1);

%% Run simulation
% simulation settings
% set compile options
configset = getconfigset(model);
compileOptions = get(configset, "CompileOptions");
set(compileOptions, "DimensionalAnalysis", false);
set(configset.SolverOptions, "MaxStep", 10)
set(configset.SolverOptions, "AbsoluteTolerance", 1e-9)
set(configset.SolverOptions, "RelativeTolerance", 1e-6)
% set simulation time
tf = 80*365*24; % 80 years in hours
set(configset, 'StopTime', tf)

% turn on estrogen option
doEST = sbioselect(model, 'Type', 'parameter', 'Name', 'do_EST');
set(doEST, 'Value', 1)


% Set age and efficacy of drug
age_ARB = 60; % years
age_ACEi = 60; % years
pctARB = 0.9359; % average value
pctACEi = 0.956; % average value

ageARB = sbioselect(model, 'Type', 'parameter', 'Name', 'age_ARB');
set(ageARB, 'Value',age_ARB)
ageACEi = sbioselect(model, 'Type', 'parameter', 'Name', 'age_ACEi');
set(ageACEi, 'Value',age_ACEi)
ARB_pct_inhib = sbioselect(model, 'Type', 'parameter', 'Name', 'pct_ARB_inhib');
set(ARB_pct_inhib, 'Value', pctARB)
ACEi_pct_inhib = sbioselect(model, 'Type', 'parameter', 'Name', 'pct_ACEi_inhib');
set(ACEi_pct_inhib, 'Value', pctACEi)

% Run simulations
fprintf('sim 1 \n')
tic
% simulation 1: estrogen decline, no ARB/ACEi
simDat1 = sbiosimulate(model);

fprintf('sim 2 \n')
% sim 2: estrogen decline + ARB at 60 years
varARB = addvariant(model, 'ARBsim');
addcontent(varARB, {'parameter', 'do_ARB', 'Value', 1});
simDat2 = sbiosimulate(model, varARB);


fprintf('sim 3 \n')
% sim 3: estrogen decline + ACEi at 60 years
varACEi = addvariant(model, 'ACEisim');
addcontent(varACEi, {'parameter', 'do_ACEi', 'Value', 1});
simDat3 = sbiosimulate(model, varACEi);

toc

%% Make figures
fprintf("plotting results \n")
% fig specs
lw = 5;
cmap = parula(8);
c1 = cmap(1,:);
c2 = cmap(4,:);
c3 = cmap(7,:);

ls_noRASi = '-';
ls_ARB = ':';
ls_ACEi = '-.';

xminmax = [40,80];
t_yrs1 = (simDat1.Time)/(365*24); % time in years
t_yrs2 = (simDat2.Time)/(365*24);
t_yrs3 = (simDat3.Time)/(365*24);

xlab = 't (years)';
fsize = 16;
labs = {'no RASi', 'ARB', 'ACEi'};
alph1 = 0.5;
alph2 = 0.6;
ms = 15;

% Bone figure
f = figure(1);
clf;
width = (2/3)*1600;
height = 900;
f.Position = [100, 100, width, height];
tiledlayout(2,2)

% BMD
nexttile(1);
id = 36;
hold on;
yline(100,"linewidth",2,"HandleVisibility", "off")
plot(t_yrs1, simDat1.Data(:,id)*100, 'linewidth',lw,...
                                    'linestyle',ls_noRASi,...
                                    'color',c1)
plot(t_yrs2, simDat2.Data(:,id)*100, 'linewidth',lw,...
                                    'linestyle',ls_ARB,...
                                    'color',c2)
plot(t_yrs3, simDat3.Data(:,id)*100, 'linewidth',lw,...
                                    'linestyle',ls_ACEi,...
                                    'color',c3)
grid on
xlim(xminmax)
xlabel(xlab)
ylabel('BMD (%)')
set(gca,'fontsize',fsize)
legend(labs,'location','southwest')
ylim([50,120])

% OC (normalized)
nexttile(2);
id = 21;
hold on
plot(t_yrs1, simDat1.Data(:,id)./model.species(id).Value * 100,...
                    'linewidth',lw,...
                    'linestyle',ls_noRASi,...
                    'color',c1)
plot(t_yrs2, simDat2.Data(:,id)./model.species(id).Value * 100,...
                    'linewidth',lw,...
                    'linestyle',ls_ARB,...
                    'color',c2)
plot(t_yrs3, simDat3.Data(:,id)./model.species(id).Value * 100,...
                    'linewidth',lw,...
                    'linestyle',ls_ACEi,...
                    'color',c3)
yline(100,"linewidth",2,"HandleVisibility","off")
grid on
xlim(xminmax)
xlabel(xlab)
set(gca,"fontsize",fsize)
legend(labs,"location","northwest")
ylabel("Osteoclasts (%)")
ylim([0,600])

% RANK-RANKL (normalized)
nexttile(3);
id = 26;
hold on
plot(t_yrs1, simDat1.Data(:,id)./model.species(id).Value * 100,...
                    'linewidth',lw,...
                    'linestyle',ls_noRASi,...
                    'color',c1)
plot(t_yrs2, simDat2.Data(:,id)./model.species(id).Value * 100,...
                    'linewidth',lw,...
                    'linestyle',ls_ARB,...
                    'color',c2)
plot(t_yrs3, simDat3.Data(:,id)./model.species(id).Value * 100,...
                    'linewidth',lw,...
                    'linestyle',ls_ACEi,...
                    'color',c3)
yline(100,"linewidth",2,"HandleVisibility","off")
grid on
xlim(xminmax)
xlabel(xlab)
set(gca,"fontsize",fsize)
legend(labs,"location","northwest")
ylabel("RANK-RANKL (%)")
ylim([0,300])

% RANK-OPG (normalized)
nexttile(4);
id = 32;
hold on
plot(t_yrs1, simDat1.Data(:,id)./model.species(id).Value * 100,...
                    'linewidth',lw,...
                    'linestyle',ls_noRASi,...
                    'color',c1)
plot(t_yrs2, simDat2.Data(:,id)./model.species(id).Value * 100,...
                    'linewidth',lw,...
                    'linestyle',ls_ARB,...
                    'color',c2)
plot(t_yrs3, simDat3.Data(:,id)./model.species(id).Value * 100,...
                    'linewidth',lw,...
                    'linestyle',ls_ACEi,...
                    'color',c3)
yline(100,"linewidth",2,"HandleVisibility","off")
grid on
xlim(xminmax)
xlabel(xlab)
set(gca,"fontsize",fsize)
legend(labs,"location","northwest")
ylabel("RANK-OPG (%)")
ylim([0,700])

% RAS model
f = figure(2);
clf;
width = 1600;
height = 900;
f.Position = [100, 100, width, height];
tiledlayout(2,3);

% AGT
nexttile(1);
id = 8;
hold on
plot(t_yrs1, simDat1.Data(:,id)/1000,...
                    'linewidth',lw,...
                    'linestyle',ls_noRASi,...
                    'color',c1)
plot(t_yrs2, simDat2.Data(:,id)/1000,...
                    'linewidth',lw,...
                    'linestyle',ls_ARB,...
                    'color',c2)
plot(t_yrs3, simDat3.Data(:,id)/1000,...
                    'linewidth',lw,...
                    'linestyle',ls_ACEi,...
                    'color',c3)
grid on
xlim(xminmax)
xlabel(xlab)
set(gca,"fontsize",fsize)
legend(labs,"location","southwest")
ylabel("[AGT] (fmol/L)")
ylim([100,600])


% Renin
nexttile(2);
id = 7;
hold on
plot(t_yrs1, simDat1.Data(:,id),...
                    'linewidth',lw,...
                    'linestyle',ls_noRASi,...
                    'color',c1)
plot(t_yrs2, simDat2.Data(:,id),...
                    'linewidth',lw,...
                    'linestyle',ls_ARB,...
                    'color',c2)
plot(t_yrs3, simDat3.Data(:,id),...
                    'linewidth',lw,...
                    'linestyle',ls_ACEi,...
                    'color',c3)
grid on
xlim(xminmax)
xlabel(xlab)
set(gca,"fontsize",fsize)
legend(labs,"location","northwest")
ylabel("[PRC] (mU/L)")
ylim([0,300])

% Ang I
nexttile(3);
id = 9;
hold on
plot(t_yrs1, simDat1.Data(:,id),...
                    'linewidth',lw,...
                    'linestyle',ls_noRASi,...
                    'color',c1)
plot(t_yrs2, simDat2.Data(:,id),...
                    'linewidth',lw,...
                    'linestyle',ls_ARB,...
                    'color',c2)
plot(t_yrs3, simDat3.Data(:,id),...
                    'linewidth',lw,...
                    'linestyle',ls_ACEi,...
                    'color',c3)
grid on
xlim(xminmax)
xlabel(xlab)
set(gca,"fontsize",fsize)
legend(labs,"location","northwest")
ylabel("[Ang I] (pmol/L)")
ylim([0,60])

% Ang II
nexttile(4);
id = 10;
hold on
plot(t_yrs1, simDat1.Data(:,id),...
                    'linewidth',lw,...
                    'linestyle',ls_noRASi,...
                    'color',c1)
plot(t_yrs2, simDat2.Data(:,id),...
                    'linewidth',lw,...
                    'linestyle',ls_ARB,...
                    'color',c2)
plot(t_yrs3, simDat3.Data(:,id),...
                    'linewidth',lw,...
                    'linestyle',ls_ACEi,...
                    'color',c3)
grid on
xlim(xminmax)
xlabel(xlab)
set(gca,"fontsize",fsize)
legend(labs,"location","northwest")
ylabel("[Ang II] (pmol/L)")
ylim([0,60])

% AT1R
nexttile(5);
id = 13;
hold on
plot(t_yrs1, simDat1.Data(:,id),...
                    'linewidth',lw,...
                    'linestyle',ls_noRASi,...
                    'color',c1)
plot(t_yrs2, simDat2.Data(:,id),...
                    'linewidth',lw,...
                    'linestyle',ls_ARB,...
                    'color',c2)
plot(t_yrs3, simDat3.Data(:,id),...
                    'linewidth',lw,...
                    'linestyle',ls_ACEi,...
                    'color',c3)
grid on
xlim(xminmax)
xlabel(xlab)
set(gca,"fontsize",fsize)
legend(labs,"location","northwest")
ylabel("[AT1R-bound Ang II] (pmol/L)")
ylim([0,15])

% AT1R
nexttile(6);
id = 14;
hold on
plot(t_yrs1, simDat1.Data(:,id),...
                    'linewidth',lw,...
                    'linestyle',ls_noRASi,...
                    'color',c1)
plot(t_yrs2, simDat2.Data(:,id),...
                    'linewidth',lw,...
                    'linestyle',ls_ARB,...
                    'color',c2)
plot(t_yrs3, simDat3.Data(:,id),...
                    'linewidth',lw,...
                    'linestyle',ls_ACEi,...
                    'color',c3)
grid on
xlim(xminmax)
xlabel(xlab)
set(gca,"fontsize",fsize)
legend(labs,"location","northwest")
ylabel("[AT2R-bound Ang II] (pmol/L)")
ylim([0,30])


