% Simulations of estrogen decline with RASi

clear all;

tic
tvals = 20:0.25:80; % time values to output

% Run simulations
% sim 1: no RASi
fprintf('no RASi sims \n')
y_noRASi = do_RASi_sim(1, 1, 0, 0, 0, 0, 1, tvals);

% sims 2: ARB
fprintf('ARB sims \n')
pct_ARB_vals = sort([92,94.7,93.7,89,96,88,97,98.3])./100; % pct inhibition based on Hallow et al 2014
age_ARB = 60;

Y_ARB = cell(size(pct_ARB_vals));
for ii = 1:length(pct_ARB_vals)
    pct_ARB = pct_ARB_vals(ii);
    %disp(pct_ARB)
    y = do_RASi_sim(1, 0, pct_ARB, 0, age_ARB, 0, 1,tvals);
    Y_ARB{ii} = y;
end

% sims 3: ACEi
fprintf('ACEi sims \n')
pct_ACEi_vals = sort([97.2, 94])./100; % pct inhibition based on Hallow et al 2014
age_ACEi = 60;

Y_ACEi = cell(size(pct_ACEi_vals));
for ii = 1:length(pct_ACEi_vals)
    pct_ACEi = pct_ACEi_vals(ii);
    %disp(pct_ACEi)
    y = do_RASi_sim(0, 1, 0, pct_ACEi, 0, age_ACEi, 1, tvals);
    Y_ACEi{ii} = y;
end
toc

%% plot results
fprintf('plotting results \n')
% fig specs
lw = 5;
cmap = parula(12);

ls_noRASi = '-';
ls_ARB = ':';
ls_ACEi = '-.';

cmap2 = spring(3);

cmap3 = parula(8);
c_noRASi = cmap3(1,:);
c_ARB = cmap3(4,:);
c_ACEi = cmap3(7,:);

xminmax=[40,80];
xlab = 't (years)';
fsize = 16;
labs = {'no RASi', 'ARB', 'ACEi'};
alph1 = 0.5;
alph2 = 0.6;
ms = 15;

% test individual points
% figure(1)
clf;
id = 32; % at1r
hold on
plot(tvals, y_noRASi(:,id), 'linewidth', lw,...
                                        'color', c_noRASi)
for ii = 1:length(pct_ARB_vals)
    y = Y_ARB{ii};
    
    c = cmap(ii + 3,:);
    if or(ii == 1, ii == 8)
        ls = '--';
    else
        ls = '-';
    end
    plot(tvals,y(:,id),'linewidth',lw,'color',c, 'linestyle', ls)
end

for ii = 1:length(pct_ACEi_vals)
    y = Y_ACEi{ii};

    c = cmap2(ii, :);
    plot(tvals,y(:,id), 'linewidth', lw,'color',c)
end

% RAS results
f = figure(2);
clf;
width = 1600;
height = 900;
f.Position = [100, 100, width, height];
tiledlayout(2,3);

% AGT
nexttile(1);
id = 27;
hold on
plot(tvals, y_noRASi(:,id)/1000, 'linewidth', lw,...
                                'linestyle',ls_noRASi,...
                                'color', c_noRASi)
[minARB, maxARB, meanARB] = get_vals(Y_ARB, id, tvals);
fill([tvals, fliplr(tvals)], [minARB./1000, fliplr(maxARB./1000)],...
                        c_ARB,...
                        'FaceColor', c_ARB,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none','HandleVisibility', 'off')
plot(tvals, meanARB/1000, 'linewidth', lw, 'color',c_ARB,...
                    'linestyle', ls_ARB)
[minACEi, maxACEi, meanACEi] = get_vals(Y_ACEi,id,tvals);
fill([tvals, fliplr(tvals)], [minACEi./1000, fliplr(maxACEi./1000)],...
                        c_ACEi,...
                        'FaceColor', c_ACEi,...
                        'FaceAlpha', alph2, ...
                        'EdgeColor', 'none',...
                        'HandleVisibility', 'off')
plot(tvals, meanACEi/1000, 'linewidth', lw, 'color',c_ACEi,...
                        'linestyle',ls_ACEi)
grid on
xlim(xminmax)
ylim([100,600])
xlabel(xlab)
ylabel('[AGT] (fmol/L)')
set(gca,'fontsize',fsize)
legend(labs,'location','southwest')

% Renin (PRC)
nexttile(2);
id = 26;
hold on
plot(tvals, y_noRASi(:,id), 'linewidth', lw,...
                                'linestyle',ls_noRASi,...
                                'color', c_noRASi)
[minARB, maxARB, meanARB] = get_vals(Y_ARB, id, tvals);
fill([tvals, fliplr(tvals)], [minARB, fliplr(maxARB)],...
                        c_ARB,...
                        'FaceColor', c_ARB,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none','HandleVisibility', 'off')
plot(tvals, meanARB, 'linewidth', lw, 'color',c_ARB, 'linestyle', ls_ARB)
[minACEi, maxACEi, meanACEi] = get_vals(Y_ACEi,id,tvals);
fill([tvals, fliplr(tvals)], [minACEi, fliplr(maxACEi)],...
                        c_ACEi,...
                        'FaceColor', c_ACEi,...
                        'FaceAlpha', alph2, ...
                        'EdgeColor', 'none',...
                        'HandleVisibility', 'off')
plot(tvals, meanACEi, 'linewidth', lw, 'color',c_ACEi, 'linestyle', ls_ACEi)
grid on
xlim(xminmax)
ylim([0,325])
xlabel(xlab)
ylabel('[Renin] (mU/L)')
set(gca,'fontsize',fsize)
legend(labs,'location','northwest')

% Ang I
nexttile(3);
id = 28;
hold on
plot(tvals, y_noRASi(:,id), 'linewidth', lw,...
                                'linestyle',ls_noRASi,...
                                'color', c_noRASi)
[minARB, maxARB, meanARB] = get_vals(Y_ARB, id, tvals);
fill([tvals, fliplr(tvals)], [minARB, fliplr(maxARB)],...
                        c_ARB,...
                        'FaceColor', c_ARB,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none','HandleVisibility', 'off')
plot(tvals, meanARB, 'linewidth', lw, 'color',c_ARB, 'linestyle', ls_ARB)
[minACEi, maxACEi, meanACEi] = get_vals(Y_ACEi,id,tvals);
fill([tvals, fliplr(tvals)], [minACEi, fliplr(maxACEi)],...
                        c_ACEi,...
                        'FaceColor', c_ACEi,...
                        'FaceAlpha', alph2, ...
                        'EdgeColor', 'none',...
                        'HandleVisibility', 'off')
plot(tvals, meanACEi, 'linewidth', lw, 'color',c_ACEi, 'linestyle', ls_ACEi)
grid on
xlim(xminmax)
ylim([0,40])
xlabel(xlab)
ylabel('[Ang I] (pmol/L)')
set(gca,'fontsize',fsize)
legend(labs,'location','northwest')

% Ang II
nexttile(4);
id = 29;
hold on
plot(tvals, y_noRASi(:,id), 'linewidth', lw,...
                                'linestyle',ls_noRASi,...
                                'color', c_noRASi)
[minARB, maxARB, meanARB] = get_vals(Y_ARB, id, tvals);
fill([tvals, fliplr(tvals)], [minARB, fliplr(maxARB)],...
                        c_ARB,...
                        'FaceColor', c_ARB,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none','HandleVisibility', 'off')
plot(tvals, meanARB, 'linewidth', lw, 'color',c_ARB, 'linestyle', ls_ARB)
[minACEi, maxACEi, meanACEi] = get_vals(Y_ACEi,id,tvals);
fill([tvals, fliplr(tvals)], [minACEi, fliplr(maxACEi)],...
                        c_ACEi,...
                        'FaceColor', c_ACEi,...
                        'FaceAlpha', alph2, ...
                        'EdgeColor', 'none',...
                        'HandleVisibility', 'off')
plot(tvals, meanACEi, 'linewidth', lw, 'color',c_ACEi, 'linestyle', ls_ACEi)
grid on
xlim(xminmax)
ylim([0,70])
xlabel(xlab)
ylabel('[Ang II] (pmol/L)')
set(gca,'fontsize',fsize)
legend(labs,'location','northwest')

% AT1R
nexttile(5);
id = 32;
hold on
plot(tvals, y_noRASi(:,id), 'linewidth', lw,...
                                'linestyle',ls_noRASi,...
                                'color', c_noRASi)
[minARB, maxARB, meanARB] = get_vals(Y_ARB, id, tvals);
fill([tvals, fliplr(tvals)], [minARB, fliplr(maxARB)],...
                        c_ARB,...
                        'FaceColor', c_ARB,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none','HandleVisibility', 'off')
plot(tvals, meanARB, 'linewidth', lw, 'color',c_ARB, 'linestyle', ls_ARB)
[minACEi, maxACEi, meanACEi] = get_vals(Y_ACEi,id,tvals);
fill([tvals, fliplr(tvals)], [minACEi, fliplr(maxACEi)],...
                        c_ACEi,...
                        'FaceColor', c_ACEi,...
                        'FaceAlpha', alph2, ...
                        'EdgeColor', 'none',...
                        'HandleVisibility', 'off')
plot(tvals, meanACEi, 'linewidth', lw, 'color',c_ACEi, 'linestyle', ls_ACEi)
grid on
xlim(xminmax)
ylim([0,15])
xlabel(xlab)
ylabel('[AT1R-bound Ang II] (pmol/L)')
set(gca,'fontsize',fsize)
legend(labs,'location','northwest')

% AT2R
nexttile(6);
id = 33;
hold on
plot(tvals, y_noRASi(:,id), 'linewidth', lw,...
                                'linestyle',ls_noRASi,...
                                'color', c_noRASi)
[minARB, maxARB, meanARB] = get_vals(Y_ARB, id, tvals);
fill([tvals, fliplr(tvals)], [minARB, fliplr(maxARB)],...
                        c_ARB,...
                        'FaceColor', c_ARB,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none','HandleVisibility', 'off')
plot(tvals, meanARB, 'linewidth', lw, 'color',c_ARB, 'linestyle', ls_ARB)
[minACEi, maxACEi, meanACEi] = get_vals(Y_ACEi,id,tvals);
fill([tvals, fliplr(tvals)], [minACEi, fliplr(maxACEi)],...
                        c_ACEi,...
                        'FaceColor', c_ACEi,...
                        'FaceAlpha', alph2, ...
                        'EdgeColor', 'none',...
                        'HandleVisibility', 'off')
plot(tvals, meanACEi, 'linewidth', lw, 'color',c_ACEi, 'linestyle', ls_ACEi)
grid on
xlim(xminmax)
ylim([0,40])
xlabel(xlab)
ylabel('[AT2R-bound Ang II] (pmol/L)')
set(gca,'fontsize',fsize)
legend(labs,'location','northwest')

AddLetters2Plots(figure(2), {'(A)', '(B)', '(C)', '(D)', '(E)','(F)'},...
                'HShift', -0.06, 'VShift', -0.06,...
                'fontsize', 20)

% Bone figure
f = figure(3);
clf; 
width = (2/3)*1600;
height = 900;
f.Position = [100,100,width,height];
tiledlayout(2,2);

% BMD
nexttile(1);
id = 34;
hold on
plot(tvals, y_noRASi(:,id)*100, 'linewidth', lw,...
                                'linestyle',ls_noRASi,...
                                'color', c_noRASi)
% load looker 1998 data
T_BMD = load("./data/Looker1998_BMDfn_female.mat").T;
errorbar(T_BMD.Age,...
    100*T_BMD.Mean_all./T_BMD.Mean_all(1), ...
    100*T_BMD.SD_all./T_BMD.Mean_all(1),...
                    'marker', 'o',...
                    'markerfacecolor',c_noRASi,...
                    'markersize', ms,...
                    'color', c_noRASi,...
                    'linewidth',2,'linestyle','none',...
                    'handlevisibility','off')
[minARB, maxARB, meanARB] = get_vals(Y_ARB, id, tvals);
fill([tvals, fliplr(tvals)], [minARB*100, fliplr(maxARB*100)],...
                        c_ARB,...
                        'FaceColor', c_ARB,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none','HandleVisibility', 'off')
plot(tvals, meanARB*100, 'linewidth', lw, 'color',c_ARB, 'linestyle', ls_ARB)
[minACEi, maxACEi, meanACEi] = get_vals(Y_ACEi,id,tvals);
fill([tvals, fliplr(tvals)], [minACEi*100, fliplr(maxACEi*100)],...
                        c_ACEi,...
                        'FaceColor', c_ACEi,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none',...
                        'HandleVisibility', 'off')
plot(tvals, meanACEi*100, 'linewidth', lw, 'color',c_ACEi, 'linestyle', ls_ACEi)
grid on
xlim([20,80])
ylim([50,120])
xlabel(xlab)
ylabel('Bone mineral density (%)')
set(gca,'fontsize',fsize)
legend(labs,'location','southwest')

% OC (normalized) 
IC = y_noRASi(1,:);
nexttile(2);
id = 11;
hold on
plot(tvals, y_noRASi(:,id)./IC(id)*100, 'linewidth', lw,...
                                'linestyle',ls_noRASi,...
                                'color', c_noRASi)
[minARB, maxARB, meanARB] = get_vals(Y_ARB, id, tvals);
% normalize
minARB = minARB./IC(id)*100;
maxARB = maxARB./IC(id)*100;
meanARB = meanARB./IC(id)*100;
fill([tvals, fliplr(tvals)], [minARB, fliplr(maxARB)],...
                        c_ARB,...
                        'FaceColor', c_ARB,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none','HandleVisibility', 'off')
plot(tvals, meanARB, 'linewidth', lw, 'color',c_ARB, 'linestyle', ls_ARB)
[minACEi, maxACEi, meanACEi] = get_vals(Y_ACEi,id,tvals);
% normalize
minACEi = minACEi./IC(id)*100;
maxACEi = maxACEi./IC(id)*100;
meanACEi = meanACEi./IC(id)*100;
fill([tvals, fliplr(tvals)], [minACEi, fliplr(maxACEi)],...
                        c_ACEi,...
                        'FaceColor', c_ACEi,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none',...
                        'HandleVisibility', 'off')
plot(tvals, meanACEi, 'linewidth', lw, 'color',c_ACEi, 'linestyle', ls_ACEi)
grid on
xlabel(xlab)
ylim([0,600])
xlim(xminmax)
ylabel('Osteoclasts (%)')
set(gca,'fontsize',fsize)
legend(labs,'location','northwest')

% RANK-RANKL (normalized) 
nexttile(3);
id = 14;
hold on
plot(tvals, y_noRASi(:,id)./IC(id)*100, 'linewidth', lw,...
                                'linestyle',ls_noRASi,...
                                'color', c_noRASi)
[minARB, maxARB, meanARB] = get_vals(Y_ARB, id, tvals);
% normalize
minARB = minARB./IC(id)*100;
maxARB = maxARB./IC(id)*100;
meanARB = meanARB./IC(id)*100;
fill([tvals, fliplr(tvals)], [minARB, fliplr(maxARB)],...
                        c_ARB,...
                        'FaceColor', c_ARB,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none','HandleVisibility', 'off')
plot(tvals, meanARB, 'linewidth', lw, 'color',c_ARB, 'linestyle', ls_ARB)
[minACEi, maxACEi, meanACEi] = get_vals(Y_ACEi,id,tvals);
% normalize
minACEi = minACEi./IC(id)*100;
maxACEi = maxACEi./IC(id)*100;
meanACEi = meanACEi./IC(id)*100;
fill([tvals, fliplr(tvals)], [minACEi, fliplr(maxACEi)],...
                        c_ACEi,...
                        'FaceColor', c_ACEi,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none',...
                        'HandleVisibility', 'off')
plot(tvals, meanACEi, 'linewidth', lw, 'color',c_ACEi, 'linestyle', ls_ACEi)
grid on
xlabel(xlab)
xlim(xminmax)
ylim([0,300])
ylabel('RANK-RANKL (%)')
set(gca,'fontsize',fsize)
legend(labs,'location','northwest')

% RANKL-OPG (normalized) 
nexttile(4);
id = 22;
hold on
plot(tvals, y_noRASi(:,id)./IC(id)*100, 'linewidth', lw,...
                                'linestyle',ls_noRASi,...
                                'color', c_noRASi)
[minARB, maxARB, meanARB] = get_vals(Y_ARB, id, tvals);
% normalize
minARB = minARB./IC(id)*100;
maxARB = maxARB./IC(id)*100;
meanARB = meanARB./IC(id)*100;
fill([tvals, fliplr(tvals)], [minARB, fliplr(maxARB)],...
                        c_ARB,...
                        'FaceColor', c_ARB,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none','HandleVisibility', 'off')
plot(tvals, meanARB, 'linewidth', lw, 'color',c_ARB, 'linestyle', ls_ARB)
[minACEi, maxACEi, meanACEi] = get_vals(Y_ACEi,id,tvals);
% normalize
minACEi = minACEi./IC(id)*100;
maxACEi = maxACEi./IC(id)*100;
meanACEi = meanACEi./IC(id)*100;
fill([tvals, fliplr(tvals)], [minACEi, fliplr(maxACEi)],...
                        c_ACEi,...
                        'FaceColor', c_ACEi,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none',...
                        'HandleVisibility', 'off')
plot(tvals, meanACEi, 'linewidth', lw, 'color',c_ACEi, 'linestyle', ls_ACEi)
grid on
xlabel(xlab)
xlim(xminmax)
%ylim([0,700])
ylabel('RANKL-OPG (%)')
set(gca,'fontsize',fsize)
legend(labs,'location','northwest')

AddLetters2Plots(figure(3), {'(A)', '(B)', '(C)', '(D)'},...
                'HShift', -0.06, 'VShift', -0.06,...
                'fontsize', 20)




% Calcium homeostasis
f = figure(4);
clf;
width = (2/3)*1600;
height = 900;
f.Position = [100, 100, width, height];
tiledlayout(2,2);

p = set_params();

% PTH
nexttile(1);
id = 3;
hold on
plot(tvals, y_noRASi(:,id)/p.Vp, 'linewidth', lw,...
                                'linestyle',ls_noRASi,...
                                'color', c_noRASi)
[minARB, maxARB, meanARB] = get_vals(Y_ARB, id, tvals);
fill([tvals, fliplr(tvals)], [minARB/p.Vp, fliplr(maxARB/p.Vp)],...
                        c_ARB,...
                        'FaceColor', c_ARB,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none','HandleVisibility', 'off')
plot(tvals, meanARB/p.Vp, 'linewidth', lw, 'color',c_ARB, 'linestyle', ls_ARB)
[minACEi, maxACEi, meanACEi] = get_vals(Y_ACEi,id,tvals);
fill([tvals, fliplr(tvals)], [minACEi/p.Vp, fliplr(maxACEi/p.Vp)],...
                        c_ACEi,...
                        'FaceColor', c_ACEi,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none',...
                        'HandleVisibility', 'off')
plot(tvals, meanACEi/p.Vp, 'linewidth', lw, 'color',c_ACEi, 'linestyle', ls_ACEi)
yline(1.1,'linewidth',2,"handlevisibility","off",'linestyle','--')
yline(6.1,'linewidth',2,"handlevisibility","off",'linestyle','--')
grid on
xlim(xminmax)
xlabel(xlab)
ylim([1,7])
ylabel('[PTH] (pmol/L)')
set(gca,'fontsize',fsize)
legend(labs,'location','southeast')

% Calcitriol
nexttile(2);
id = 5;
hold on
plot(tvals, y_noRASi(:,id)/p.Vp, 'linewidth', lw,...
                                'linestyle',ls_noRASi,...
                                'color', c_noRASi)
[minARB, maxARB, meanARB] = get_vals(Y_ARB, id, tvals);
fill([tvals, fliplr(tvals)], [minARB/p.Vp, fliplr(maxARB/p.Vp)],...
                        c_ARB,...
                        'FaceColor', c_ARB,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none','HandleVisibility', 'off')
plot(tvals, meanARB/p.Vp, 'linewidth', lw, 'color',c_ARB, 'linestyle', ls_ARB)
[minACEi, maxACEi, meanACEi] = get_vals(Y_ACEi,id,tvals);
fill([tvals, fliplr(tvals)], [minACEi/p.Vp, fliplr(maxACEi/p.Vp)],...
                        c_ACEi,...
                        'FaceColor', c_ACEi,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none',...
                        'HandleVisibility', 'off')
plot(tvals, meanACEi/p.Vp, 'linewidth', lw, 'color',c_ACEi, 'linestyle', ls_ACEi)
yline(60,'linewidth',2,"handlevisibility","off",'linestyle','--')
yline(195,'linewidth',2,"handlevisibility","off",'linestyle','--')
grid on
xlim(xminmax)
xlabel(xlab)
ylim([50,200])
ylabel('[Calcitriol] (pmol/L)')
set(gca,'fontsize',fsize)
legend(labs,'location','best')

% Calcium
nexttile(3);
id = 7;
hold on
plot(tvals, y_noRASi(:,id)/p.Vp, 'linewidth', lw,...
                                'linestyle',ls_noRASi,...
                                'color', c_noRASi)
[minARB, maxARB, meanARB] = get_vals(Y_ARB, id, tvals);
fill([tvals, fliplr(tvals)], [minARB/p.Vp, fliplr(maxARB/p.Vp)],...
                        c_ARB,...
                        'FaceColor', c_ARB,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none','HandleVisibility', 'off')
plot(tvals, meanARB/p.Vp, 'linewidth', lw, 'color',c_ARB, 'linestyle', ls_ARB)
[minACEi, maxACEi, meanACEi] = get_vals(Y_ACEi,id,tvals);
fill([tvals, fliplr(tvals)], [minACEi/p.Vp, fliplr(maxACEi/p.Vp)],...
                        c_ACEi,...
                        'FaceColor', c_ACEi,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none',...
                        'HandleVisibility', 'off')
plot(tvals, meanACEi/p.Vp, 'linewidth', lw, 'color',c_ACEi, 'linestyle', ls_ACEi)
yline(2.2,'linewidth',2,"handlevisibility","off",'linestyle','--')
yline(2.6,'linewidth',2,"handlevisibility","off",'linestyle','--')
grid on
xlim(xminmax)
xlabel(xlab)
ylim([2,3])
ylabel('[Calcium] (mmol/L)')
set(gca,'fontsize',fsize)
legend(labs,'location','northwest')

% Phosphte
nexttile(4);
id = 8;
hold on
plot(tvals, y_noRASi(:,id)/p.Vp, 'linewidth', lw,...
                                'linestyle',ls_noRASi,...
                                'color', c_noRASi)
[minARB, maxARB, meanARB] = get_vals(Y_ARB, id, tvals);
fill([tvals, fliplr(tvals)], [minARB/p.Vp, fliplr(maxARB/p.Vp)],...
                        c_ARB,...
                        'FaceColor', c_ARB,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none','HandleVisibility', 'off')
plot(tvals, meanARB/p.Vp, 'linewidth', lw, 'color',c_ARB, 'linestyle', ls_ARB)
[minACEi, maxACEi, meanACEi] = get_vals(Y_ACEi,id,tvals);
fill([tvals, fliplr(tvals)], [minACEi/p.Vp, fliplr(maxACEi/p.Vp)],...
                        c_ACEi,...
                        'FaceColor', c_ACEi,...
                        'FaceAlpha', alph1, ...
                        'EdgeColor', 'none',...
                        'HandleVisibility', 'off')
plot(tvals, meanACEi/p.Vp, 'linewidth', lw, 'color',c_ACEi, 'linestyle', ls_ACEi)
yline(0.9,'linewidth',2,"handlevisibility","off",'linestyle','--')
yline(1.5,'linewidth',2,"handlevisibility","off",'linestyle','--')
grid on
xlim(xminmax)
ylim([0.5,2])
xlabel(xlab)
ylabel('[Phosphate] (mmol/L)')
set(gca,'fontsize',fsize)
legend(labs,'location','northwest')

AddLetters2Plots(figure(4), {'(A)', '(B)', '(C)','(D)'},...
                'HShift', -0.06, 'VShift', -0.06,...
                'fontsize', 20)


%%% functions used
function [minvals,maxvals,meanvals] = get_vals(Y, id, tvals)
vals = nan(length(Y), length(tvals));
for ii = 1:length(Y)
    y = Y{ii};
    vals(ii, :) = y(:,id);
end
minvals = min(vals,[],1);
maxvals = max(vals,[],1);
meanvals = mean(vals,1);
end
