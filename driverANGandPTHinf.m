% Script for running ANG II and PTH infusion experiments
clearvars;

% TERI dose
TERI_DOSEmcg = 50; % dose in mcg
TERI_DOSE    = TERI_DOSEmcg*1e6/4117.8; % convert dose to pmol

% Set parameters
p = set_params();
% change pars here


[params, parnames] = pars2vector(p,0);

% set initial conditions
IC = load("IC/2025-04-18_ICfinal3.mat").IC;

% set time span
t0 = -2;
tf = 40; 

%% ANG II inf only
do_ANGII_inf = 1;
do_TERIiv = 0; % TERI iv dose

% Run simulation
fprintf('sim 1\n')

tspan = [t0,tf];
opts_ode = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9, 'MaxStep', 0.1);
[t1, y1] = ode15s(@(t,y) mod_eqns(t,y,params,...
                            'do_ANGII_inf', do_ANGII_inf),...
                            tspan, IC, opts_ode);


%% PTH inf only (TERI SC)
fprintf('sim 2 \n')
% pre-dosing
tspan = [t0, 0]; % pre-dosing

[t2_0,y2_0] = ode15s(@(t,y) mod_eqns(t,y,params),...
                            tspan, IC, opts_ode);


IC1 = y2_0(end,:);
IC1(10) = TERI_DOSE;

tspan = [0,tf];
[t2_1,y2_1] = ode15s(@(t,y) mod_eqns(t,y,params),...
                            tspan, IC1, opts_ode);
t2 = [t2_0;t2_1];
y2 = [y2_0;y2_1];


%% Ang II + PTH inf
do_ANGII_inf = 1;

% Run simulation
fprintf('sim 3 \n')


tspan = [t0,0];
opts_ode = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9, 'MaxStep', 0.1);
[t3_0, y3_0] = ode15s(@(t,y) mod_eqns(t,y,params,...
                            'do_ANGII_inf', do_ANGII_inf),...
                            tspan, IC, opts_ode);
IC1 = y3_0(end,:);
IC1(10) = TERI_DOSE;

tspan = [0,tf];
[t3_1,y3_1] = ode15s(@(t,y) mod_eqns(t,y,params,...
                            'do_ANGII_inf', do_ANGII_inf),...
                            tspan, IC1, opts_ode);
t3 = [t3_0;t3_1];
y3 = [y3_0;y3_1];



%% Plot results
fprintf('plotting results. \n')

% fig specs
lw = 5;
fsize = 16;
fleg = 16;
ft = 20;
cgray = uint8([70 78 81]);
cmap = parula(10);
c1 = cmap(1,:);
c2 = cmap(6,:);
c3 = cmap(9,:);
ls1 = '-';
ls2 = '-';
ls3 = '-.';
xlab = 't (hours)';
ms = 10;
lw_err=2;
cm = cmap(5,:);
xminmax = [0,4]; %[0,24]; 
labs = {'Ang II only', 'PTH only', 'Ang II + PTH'};

% data
% load Grant 1992 data
T_PRA = load("./data/Fig3_Grant_PRA.mat").T_PRA;
T_PTH = load("./data/Fig2_Grant_PTH.mat").T_PTH;
T_AII = load("./data/AngII_Grant.mat").T_AngII;

% manuscript figure
fig = figure(1);
clf;

% Set the figure position and size
width = 1250;
height = 900;
set(fig, 'Position', [100, 100, width, height]);
% 
tiledlayout(2,2);

% PRA
nexttile(1);
id = 26;
hold on
cf = 0.1828/16.63;
plot(t1,y1(:,id)*cf,'linewidth',lw, 'color', c1, 'linestyle', ls1)
plot(t2,y2(:,id)*cf,'linewidth',lw, 'color', c2, 'linestyle', ls2)
plot(t3,y3(:,id)*cf,'linewidth',lw, 'color', c3, 'linestyle', ls3)
% plot Grant 1992 data
errorbar(T_PRA.Time/60,...
        T_PRA.AIIinf,...
        T_PRA.AIIinf_err,...
        'linestyle','none',...
        'linewidth',lw_err,...
        'marker', 'o',...
        'markersize',ms,...
        'color',c1,'MarkerFaceColor', c1,...
        'HandleVisibility', 'off')
errorbar(T_PRA.Time/60,...
        T_PRA.PTHinf,...
        T_PRA.PTHinf_err,...
        'linestyle','none',...
        'linewidth',lw_err,...
        'marker', 'o',...
        'markersize',ms,...
        'color',c2,'MarkerFaceColor', c2,...
        'HandleVisibility', 'off')
errorbar(T_PRA.Time/60,...
        T_PRA.BothInf,...
        T_PRA.Bothinf_err,...
        'linestyle','none',...
        'linewidth',lw_err,...
        'marker', 'o',...
        'markersize',ms,...
        'color',c3,'MarkerFaceColor', c3,...
        'HandleVisibility', 'off')

legend(labs, 'location', 'northeast', 'fontsize', fleg)
xlabel(xlab)
ylim([0,0.6])
xlim(xminmax)
ylabel('PRA (ng/(L\cdot s))') % PRA
set(gca,'fontsize',fsize)
grid on


% Ang II
nexttile(2);
id = 29;
hold on
plot(t1,y1(:,id),'linewidth',lw, 'color', c1, 'linestyle', ls1)

plot(t2,y2(:,id),'linewidth',lw, 'color', c2, 'linestyle', ls2)
plot(t3,y3(:,id),'linewidth',lw, 'color', c3, 'linestyle', ls3)
% plot Grant 1992 data
cf = y1(1,29)/T_AII.AIIinf(1);
errorbar(T_AII.Time/60,...
        cf*T_AII.AIIinf,...
        cf*T_AII.AIIinf_err,...
            'linestyle', 'none',...
            'linewidth', lw_err,...
            'marker', 'o',...
            'markersize', ms,...
            'color', c1, 'MarkerFaceColor', c1,...
            'HandleVisibility', 'off')
errorbar(T_AII.Time/60,...
        cf*T_AII.PTHinf,...
        cf*T_AII.PTHinf_err,...
            'linestyle', 'none',...
            'linewidth', lw_err,...
            'marker', 'o',...
            'markersize', ms,...
            'color', c2, 'MarkerFaceColor', c2,...
            'HandleVisibility', 'off')
errorbar(T_AII.Time/60,...
        cf*T_AII.Bothinf,...
        cf*T_AII.BothInf_err,...
            'linestyle', 'none',...
            'linewidth', lw_err,...
            'marker', 'o',...
            'markersize', ms,...
            'color', c3, 'MarkerFaceColor', c3,...
            'HandleVisibility', 'off')
legend(labs, 'location', 'northeast', 'fontsize',fleg)
xlabel(xlab)
ylim([0,40])
ylabel('[AngII] (pmol/L)')
xlim(xminmax)
set(gca,'fontsize',fsize)
grid on



% AT1R
nexttile(4);
id = 32;
hold on
plot(t1,y1(:,id),'linewidth',lw, 'color', c1, 'linestyle', ls1)
plot(t2,y2(:,id),'linewidth',lw, 'color', c2, 'linestyle', ls2)
plot(t3,y3(:,id),'linewidth',lw, 'color', c3, 'linestyle', ls3)
legend(labs, 'location', 'northeast', 'fontsize',fleg)
xlabel(xlab)
xlim(xminmax)
ylim([0,20])
ylabel('[AT1R-bound Ang II] (pmol/L)')
set(gca,'fontsize',fsize)
grid on


% PTH
nexttile(3);
id = 3;
hold on
plot(t1,y1(:,id)/p.Vp,'linewidth',lw, 'color', c1, 'linestyle', ls1)

plot(t2,y2(:,id)/p.Vp,'linewidth',lw, 'color', c2, 'linestyle', ls2)
plot(t3,y3(:,id)/p.Vp,'linewidth',lw, 'color', c3, 'linestyle', ls3)
% plot Grant 1992 data
cf = (y1(1,id)/p.Vp)/T_PTH.AIIinf(1);
errorbar(T_PTH.Time/60,...
            T_PTH.AIIinf*cf,...
            T_PTH.AIIinf_err*cf,...
            'linestyle','none',...
            'linewidth',lw_err,...
            'marker', 'o',...
            'markersize', ms,...
            'color', c1, 'MarkerFaceColor', c1,...
            'HandleVisibility', 'off')
errorbar(T_PTH.Time/60,...
            T_PTH.PTHinf*cf,...
            T_PTH.PTHinf_err*cf,...
            'linestyle','none',...
            'linewidth',lw_err,...
            'marker', 'o',...
            'markersize', ms,...
            'color', c2, 'MarkerFaceColor', c2,...
            'HandleVisibility', 'off')
errorbar(T_PTH.Time/60,...
            T_PTH.BothInf*cf,...
            T_PTH.Bothinf_err*cf,...
            'linestyle','none',...
            'linewidth',lw_err,...
            'marker', 'o',...
            'markersize', ms,...
            'color', c3, 'MarkerFaceColor', c3,...
            'HandleVisibility', 'off')
legend(labs, 'location', 'northeast','fontsize',fleg)
ylim([0,6])
xlabel(xlab)
xlim(xminmax)
ylabel('[iPTH] (pmol/L)')
set(gca,'fontsize',fsize)
grid on


% Add letters to the figure
AddLetters2Plots(figure(1), {'(A)', '(B)', '(C)', '(D)'},...
                'HShift', -0.09,'VShift',-0.02,... 
                'fontsize', 20)


% % TERI
% figure(2);
% clf;
% tiledlayout(1,2)
% 
% % TERISC
% nexttile(1);
% id = 10;
% hold on
% plot(t1,y1(:,id),'linewidth',lw, 'color', c1, 'linestyle', ls1)
% plot(t2,y2(:,id),'linewidth',lw, 'color', c2, 'linestyle', ls2)
% plot(t3,y3(:,id),'linewidth',lw, 'color', c3, 'linestyle', ls3)
% legend(labs, 'location', 'northeast', 'fontsize',fleg)
% xlabel(xlab)
% xlim(xminmax)
% ylabel('TERISC')
% set(gca,'fontsize',fsize)
% grid on
% 
% % TERI Central
% nexttile(2);
% id = 35;
% hold on
% plot(t1,y1(:,id),'linewidth',lw, 'color', c1, 'linestyle', ls1)
% plot(t2,y2(:,id),'linewidth',lw, 'color', c2, 'linestyle', ls2)
% plot(t3,y3(:,id),'linewidth',lw, 'color', c3, 'linestyle', ls3)
% legend(labs, 'location', 'northeast', 'fontsize',fleg)
% xlabel(xlab)
% xlim(xminmax)
% ylabel('TERI CENT')
% set(gca,'fontsize',fsize)
% grid on















