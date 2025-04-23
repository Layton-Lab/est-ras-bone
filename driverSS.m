% Script for solving for SS and solving ODEs starting at SS
clearvars;


%% Set initial conditions
load("IC/2025-04-18_ICfinal3.mat","IC")

%% Set parameters
p = set_params();
% change pars here


[params, parnames] = pars2vector(p,0);


%% Optional model inputs
do_PTH = 1; 
do_Ctriol = 1; 
do_Ca = 1; 
do_Phos = 1;
do_bone = 1; 
do_RAS = 1;
do_EST = 0;


%% Run simulation starting at IC
fprintf('solving ODEs. \n')

% set time span
t0 = 0;
tf = 1e6;

tspan = [t0,tf];
opts_ode = odeset('RelTol', 1.0e-9, 'AbsTol', 1e-12, 'MaxStep', 10);
[t, y] = ode15s(@(t,y) mod_eqns(t,y,params,...
                            'do_PTH', do_PTH, ...
                            'do_Ctriol', do_Ctriol, ...
                            'do_Ca', do_Ca,...
                            'do_Phos', do_Phos,...
                            'do_RAS', do_RAS,...
                            'do_bone', do_bone,...
                            'do_EST', do_EST),...
                            tspan, IC, opts_ode);
fprintf('plotting results. \n')

%% Plot results (trajectory)
% fig specs
lw = 7;
fsize = 18;
ft = 20;
cgray = uint8([70 78 81]);
cmap = summer(6);
c1 = cmap(2,:);
c2 = cmap(4,:);
xlab = 't';
all_figs = 1; % options to show all figures


%% PTH model
if or(do_PTH, all_figs)
    figure(3)
    clf;
    nr = 2; nc = 2;
    % PTmax
    subplot(nr,nc,1)
    id = 1;
    plot(t,y(:,id),'linewidth',lw,'color',c1)
    xlabel(xlab)
    ylabel('PTmax')
    ymin = max(0,round(min(y(:,id))-1));
    ymax = round(max(y(:,id))+1);
    ylim([ymin,ymax])
    set(gca,'fontsize',fsize)
    grid on
    
    % PTG
    subplot(nr,nc,2)
    id = 2;
    plot(t,y(:,id),'linewidth',lw,'color',c1)
    xlabel(xlab)
    ymin = max(0,round(min(y(:,id))-1));
    ymax = round(max(y(:,id))+1);
    ylim([ymin,ymax])
    ylabel('PTG')
    set(gca,'fontsize',fsize)
    grid on

    % PTH
    subplot(nr,nc,3)
    id = 3;
    plot(t,y(:,id),'linewidth',lw,'color',c1)
    ymin = max(0, round(min(y(:,id) - 1)));
    ymax = round(max(y(:,id) + 1));
    ylim([ymin,ymax])
    ylabel('PTH')
    xlabel(xlab)
    set(gca,'fontsize',fsize)
    grid on

    % PTHConc
    subplot(nr,nc,4)
    plot(t,y(:,id)./p.Vp, 'linewidth', lw, 'color',c1)
    yline(1.1, 'color','k', 'linestyle', '--')
    yline(6.1, 'color', 'k', 'linestyle', '--')
    ylabel('[PTH]')
%     ymin = max(0,round(min(y(:,9))-1));
%     ymax = round(max(y(:,9))+1);
%     ylim([ymin,ymax])
    xlabel(xlab)
    set(gca,'fontsize',fsize)
    grid on

    sgtitle('PTH model', 'fontsize', ft)
end

%% Calcitriol model
if or(do_Ctriol,all_figs)
    figure(4)
    clf;
    nr = 1; nc = 2;
    % AOH
    subplot(nr,nc,1)
    id = 4;
    plot(t,y(:,id),'linewidth',lw,'color',c1)
    xlabel(xlab)
    ylabel('1\alpha-hydroxylase (AOH)')
    ymin = round(min(y(:,id))-1);
    ymax = round(max(y(:,id))+1);
    ylim([ymin,ymax])
    set(gca,'fontsize',fsize)
    grid on
    
    % Calcitriol
    subplot(nr,nc,2)
    id = 5;
    plot(t,y(:,id)./p.Vp,'linewidth',lw,'color',c1)
    hold on
    yline(60, 'color', 'k', 'linestyle', '--')
    yline(195, 'color', 'k', 'linestyle', '--')
    xlabel(xlab)
    ylabel('[Ctriol] (pmol/L)')
    set(gca,'fontsize',fsize)
    grid on

    sgtitle('Calcitriol model', 'fontsize', ft)
end


%% Calcium model
if or(do_Ca,all_figs)
    figure(6)
    clf;
    % Ca
    id = 7;
    plot(t,y(:,id)/p.Vp,'linewidth',lw,'color',c1)
    yline(2.6, 'color', 'k', 'linestyle', '--')
    yline(2.2, 'color', 'k', 'linestyle', '--')
    xlabel(xlab)
    ylabel('Ca')
    set(gca,'fontsize',fsize)
    grid on

    sgtitle('Ca model', 'fontsize', ft)
end

%% Phosphate model
if or(do_Phos, all_figs)
    figure(7)
    clf;
    subplot(1,2,1)
    % ECC phos
    id = 8;
    plot(t,y(:,id)/p.Vp, 'linewidth', lw, 'color', c1)
    hold on
    yline(0.9, 'color', 'k', 'linestyle', '--')
    yline(1.5, 'color', 'k', 'linestyle', '--')
    xlabel(xlab)
    ylabel('[PhosECC] (mmol/L)')
    set(gca, 'fontsize', fsize)
    grid on
    
    subplot(1,2,2)
    % Int phos
    id = 9;
    plot(t,y(:,id), 'linewidth', lw, 'color', c1)
    xlabel(xlab)
    ylabel('PhosInt')
    ymin = round(min(y(:,id)) - 1);
    ymax = round(max(y(:,id)) + 1);
    ylim([ymin, ymax])
    set(gca, 'fontsize', fsize)
    grid on


    sgtitle('Phos model', 'fontsize', ft)
    
end % do_phos

%% Plot the bone model
figure(10);
clf;
tiledlayout(4,4)

nexttile;
id = 11;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('OC')
set(gca,'fontsize',fsize)

nexttile;
id = 12;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('OBfast')
set(gca,'fontsize',fsize)

nexttile;
id = 13;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('OBslow')
set(gca,'fontsize',fsize)

nexttile;
id = 14;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('M (RANK-RANKL)')
set(gca,'fontsize',fsize)

nexttile;
id = 15;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y(:,id)-1);
ymax = max(y(:,id)+1);
ylim([ymin,ymax])
ylabel('RNK')
set(gca,'fontsize',fsize)


nexttile;
id = 16;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('L')
set(gca,'fontsize',fsize)

nexttile;
id = 17;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('TGFBact')
set(gca,'fontsize',fsize)

nexttile;
id = 18;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('TGFB (latent)')
set(gca,'fontsize',fsize)

nexttile;
id = 19;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('ROB1')
set(gca,'fontsize',fsize)

nexttile;
id = 20;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('BCL2')
set(gca,'fontsize',fsize)

nexttile;
id = 21;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('O (OPG)')
set(gca,'fontsize',fsize)

nexttile;
id = 22;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('N (RANKL-OPG)')
set(gca,'fontsize',fsize)


nexttile;
id = 23;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('RX2')
set(gca,'fontsize',fsize)

nexttile;
id = 24;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('CREB')
set(gca,'fontsize',fsize)

nexttile;
id = 25;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('HAp')
set(gca,'fontsize',fsize)

% Plot RAS model
figure(11);
clf;
tiledlayout(2,4);

nexttile;
id = 26;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y(:,id)-1);
ymax = max(y(:,id)+1);
ylim([ymin,ymax])
ylabel('PRC')
set(gca,'fontsize',fsize)
grid on

nexttile;
id = 27;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y(:,id)-1);
ymax = max(y(:,id)+1);
ylim([ymin,ymax])
ylabel('AGT')
set(gca,'fontsize',fsize)
grid on

nexttile;
id = 28;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y(:,id)-1);
ymax = max(y(:,id)+1);
ylim([ymin,ymax])
ylabel('AngI')
set(gca,'fontsize',fsize)
grid on

nexttile;
id = 29;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y(:,id)-1);
ymax = max(y(:,id)+1);
ylim([ymin,ymax])
ylabel('AngII')
set(gca,'fontsize',fsize)
grid on

nexttile;
id = 30;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y(:,id)-1);
ymax = max(y(:,id)+1);
ylim([ymin,ymax])
ylabel('Ang17')
set(gca,'fontsize',fsize)
grid on

nexttile;
id = 31;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y(:,id)-1);
ymax = max(y(:,id)+1);
ylim([ymin,ymax])
ylabel('AngIV')
set(gca,'fontsize',fsize)
grid on

nexttile;
id = 32;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y(:,id)-1);
ymax = max(y(:,id)+1);
ylim([ymin,ymax])
ylabel('AT1R')
set(gca,'fontsize',fsize)
grid on

nexttile;
id = 33;
plot(t,y(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y(:,id)-1);
ymax = max(y(:,id)+1);
ylim([ymin,ymax])
ylabel('AT2R')
set(gca,'fontsize',fsize)
grid on

%% BMD
figure(9)
clf;

id = 34;
plot(t,y(:,id),'linewidth',lw,'color',c1)
ylim([0.85,1.1])
xlabel(xlab)
ylabel('BMDfn')
set(gca,'fontsize',fsize)
grid on

%% Find SS solution with fsolve
do_SS = 1;
if do_SS
fprintf('find SS \n')
IG = y(end,:)'; %IC; % initial guess
opts_fsolve = optimoptions('fsolve','Display', 'none',...
                                'MaxFunEvals', 1e6,...
                                'MaxIter', 1e6,...
                                'FunctionTolerance', 1e-16);
[SSdat, residual, ...
    exitflag, output] = fsolve(@(y) mod_eqns(0, y, params,...
                                 'do_PTH', do_PTH,...
                                 'do_Ctriol', do_Ctriol,...
                                 'do_Ca', do_Ca,...
                                 'do_Phos', do_Phos,...
                                 'do_RAS', do_RAS,...
                                 'do_bone', do_bone,...
                                'do_EST', do_EST),...
                                IG, opts_fsolve);
[v, id] = max(abs(residual));
fprintf('maximum residual size: %0.1d, id = %i \n', v, id)

if exitflag ~= 1
    fprintf('WARNING: SS did not converge exitflag: %i \n', exitflag)
end

% set initial condition to SS

%% Print out solution
if exitflag == 1
    % Print steady state solutions
    fprintf('SS solution \n')
        % PTH
        fprintf('PTmax:    %0.2f \n', SSdat(1))
        fprintf('PTG:    %0.2f \n', SSdat(2))
        fprintf('PTH:    %0.2f \n', SSdat(3))
        fprintf('PTHConc: %0.2f \n', SSdat(3)/p.Vp)

        % Ctriol
        fprintf('AOH:    %0.2f \n', SSdat(4))
        fprintf('Ctriol: %0.2f \n', SSdat(5))
        fprintf('CtriolConc: %0.2f \n', SSdat(5)/p.Vp)

        % Calcium
        fprintf('Ca:     %0.2f \n', SSdat(7))
        fprintf('CaConc: %0.2f \n', SSdat(7)/p.Vp)
        
        % Phosphate
        fprintf('PhosECC: %0.2f \n', SSdat(8))
        fprintf('ECC Phos Conc: %0.2f \n', SSdat(8)/p.Vp)
        fprintf('PhosInt: %0.2f \n', SSdat(9))

        % RAS
        fprintf('PRC:   %0.2f \n', SSdat(26))
        fprintf('AGT:   %0.2f \n', SSdat(27))
        fprintf('AngI:  %0.2f \n', SSdat(28))
        fprintf('AngII: %0.2f \n', SSdat(29))
        fprintf('AT1R:  %0.2f \n', SSdat(32))
        fprintf('AT2R:  %0.2f \n', SSdat(33))


    fprintf('\n')


    %% Plot steady state solution (bar plot)
    
    ranges = [2.2, 2.6; % CaConc
                0.9, 1.5; % ECCPhosConc
                1.1, 6.1; % PTHConc
                60, 195; % CtriolConc
                ];
    
    rmid = (ranges(:,1) + ranges(:,2))./2;
    
    SSConc = [SSdat(7)/p.Vp;  % CaConc
                SSdat(8)/p.Vp;  % ECCPhos
                SSdat(3)/p.Vp;  % PTH
                SSdat(5)/p.Vp]; % Ctriol
    
    % figure specs
    cmap = summer(3);
    c1 = cmap(1,:);
    c2 = cmap(2,:);
    bw1 = 0.75;
    bw2 = 0.5;
    labs = {'CaConc', 'ECCPhosConc', 'PTHConc', 'CtriolConc'};
    
    figure(100);
    clf;
    tiledlayout(2,2);
    for ii = 1:4
       nexttile;
       %disp(labs{ii})
       %figure(temp);
       %clf;
       xvals = 1;
       bar(xvals, SSConc(ii), 'facecolor', c2, 'barwidth', bw1)
       hold on
       h = errorbar(xvals, rmid(ii),...
                            (ranges(ii,1) - rmid(ii)),...
                            (ranges(ii,2)- rmid(ii)),...
                            'linestyle', 'none',...
                            'color', 'black');
       h.LineWidth = 2; % thickness of error bars
       h.CapSize = 25; % width of error bars
       xticks(xvals)
       xticklabels(labs{ii})
       temp = sprintf("%s = %0.2f", labs{ii}, SSConc(ii));
       title(temp)
       set(gca, 'fontsize', 18)
       grid on
    end

end

fprintf('done SS! \n')

% Run model again after SS
run_check = 0;
if run_check
tspan = [0,80*365*24];
SSdat(34) = 1; % reset BMD
[t1, y1] = ode15s(@(t,y) mod_eqns(t,y,params,...
                            'do_PTH', do_PTH, ...
                            'do_Ctriol', do_Ctriol, ...
                            'do_Ca', do_Ca,...
                            'do_Phos', do_Phos,...
                            'do_RAS', do_RAS,...
                            'do_bone', do_bone,...
                            'do_EST', do_EST),...
                            tspan, SSdat, opts_ode);
figure(20);
clf;
id = 34;
plot(t1/(365 * 24), y1(:,id), 'color', c2, 'linewidth', lw)
xlabel('t (years)')
ylabel('BMD')


%% PTH
figure(21);
clf;
id = 3;
plot(t1/(365 * 24), y1(:,id), 'color', c2, 'linewidth', lw)
xlabel('t (years)')
ylabel('PTH')



%% Plot the bone model
figure(22);
clf;
tiledlayout(4,4)

nexttile;
id = 11;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('OC')
set(gca,'fontsize',fsize)

nexttile;
id = 12;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('OBfast')
set(gca,'fontsize',fsize)

nexttile;
id = 13;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('OBslow')
set(gca,'fontsize',fsize)

nexttile;
id = 14;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('M (RANK-RANKL)')
set(gca,'fontsize',fsize)

nexttile;
id = 15;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y1(:,id)-1);
ymax = max(y1(:,id)+1);
ylim([ymin,ymax])
ylabel('RNK')
set(gca,'fontsize',fsize)


nexttile;
id = 16;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('L')
set(gca,'fontsize',fsize)

nexttile;
id = 17;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('TGFBact')
set(gca,'fontsize',fsize)

nexttile;
id = 18;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('TGFB (latent)')
set(gca,'fontsize',fsize)

nexttile;
id = 19;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('ROB1')
set(gca,'fontsize',fsize)

nexttile;
id = 20;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('BCL2')
set(gca,'fontsize',fsize)

nexttile;
id = 21;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('O (OPG)')
set(gca,'fontsize',fsize)

nexttile;
id = 22;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('N (RANKL-OPG)')
set(gca,'fontsize',fsize)


nexttile;
id = 23;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('RX2')
set(gca,'fontsize',fsize)

nexttile;
id = 24;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('CREB')
set(gca,'fontsize',fsize)

nexttile;
id = 25;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
% ymin = min(y(:,id)-1);
% ymax = max(y(:,id)+1);
% ylim([ymin,ymax])
ylabel('HAp')
set(gca,'fontsize',fsize)

% Plot RAS model
figure(23);
clf;
tiledlayout(2,4);

nexttile;
id = 26;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y1(:,id)-1);
ymax = max(y1(:,id)+1);
ylim([ymin,ymax])
ylabel('PRC')
set(gca,'fontsize',fsize)
grid on

nexttile;
id = 27;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y1(:,id)-1);
ymax = max(y1(:,id)+1);
ylim([ymin,ymax])
ylabel('AGT')
set(gca,'fontsize',fsize)
grid on

nexttile;
id = 28;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y1(:,id)-1);
ymax = max(y1(:,id)+1);
ylim([ymin,ymax])
ylabel('AngI')
set(gca,'fontsize',fsize)
grid on

nexttile;
id = 29;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y1(:,id)-1);
ymax = max(y1(:,id)+1);
ylim([ymin,ymax])
ylabel('AngII')
set(gca,'fontsize',fsize)
grid on

%AT1R
nexttile;
id = 32;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y1(:,id)-1);
ymax = max(y1(:,id)+1);
ylim([ymin,ymax])
ylabel('AT1R')
set(gca,'fontsize',fsize)
grid on

nexttile;
id = 33;
plot(t1,y1(:,id),'linewidth',lw)
xlabel(xlab)
ymin = min(y1(:,id)-1);
ymax = max(y1(:,id)+1);
ylim([ymin,ymax])
ylabel('AT2R')
set(gca,'fontsize',fsize)
grid on
end % run check
end
