% Make figure of the Sobol results

% clear variables
clearvars;

% load results
sRes = load("./SobolResults/21-Apr-2025_runSobol_CaRAS_notes-newpars_final.mat").sobolResults;

%% Plot simulations
cmap = parula(5);
cmap2 = turbo(6);
figure(1);
clf;
h = plotData(sRes, 'showMedian', true, 'showMean', false,...
                "Alphas", 0.05,...
                "Observables", [1,2,3,4,6,9],...
                "MedianColor", cmap2(5,:),...
                "FaceColor", cmap(2,:));

% Figure specs
bw = 1;
mst1 = 'o';
mst2 = '*';

ms1 = 8;
ms2 = 10;
fsize = 15;
fleg = 12;
labs =["S_{first}", "S_{total}"];

t_yrs = sRes.Time/ (365*24);
npars=length(sRes.SobolIndices);
nobs = length(sRes.Observables);


t_ids = [3,4,5,6,7];
%t_yrs(t_ids)

% Get max, median for each of the parameters for Sf and St
% then can sort by max and median
maxST = nan(npars, nobs); % Columns are the different observables
maxSF = nan(npars, nobs);
medST = nan(npars, nobs);
medSF = nan(npars, nobs);
meanSF = nan(npars, nobs);
meanST = nan(npars, nobs);
sdSF = nan(npars,nobs);
sdST = nan(npars,nobs);

for ii = 1:npars
    names{ii} = postprocess_name(sRes.SobolIndices(ii, 1).Parameter);
    for jj = 1:nobs % observables
        SI = sRes.SobolIndices(ii,jj);
        SFirst = SI.FirstOrder;
        STotal = SI.TotalOrder;
        maxST(ii,jj) = max(STotal);
        maxSF(ii,jj) = max(SFirst);
        medST(ii,jj)= median(STotal);
        medSF(ii,jj) = median(SFirst);
        meanST(ii,jj)= mean(STotal);
        meanSF(ii,jj) = mean(SFirst);
        sdST(ii,jj)= std(STotal);
        sdSF(ii,jj) = std(SFirst);
    end
end

% % Determine which parameters have ST > threshold for 
% observations included
th = 0.05;
obs_ids = [1,2,3,4,6,9];
ids = find(any(maxSF(:,obs_ids) > th, 2));

%length(ids)

% colors
cmap = summer(5);
c1 = cmap(1,:);
c2 = cmap(4,:);


%% 
% % Determine which parameters have ST > threshold for 
% observations included
th = 0.04;
obs_ids = [1,2];
ids = find(any(maxSF(:,obs_ids) > th, 2));

f = figure(3);
width = 1200;
height = 900;
f.Position = [100,100,width,height];
ymax = 0.5;
clf;

tiledlayout(2,1);

% BMD
nexttile(1);
obs_id = 2;
% sort values by ST
[~, ids_st] = sort(meanST(ids, obs_id), 'descend');

ids_sort = ids(ids_st);

dat = [meanSF(ids_sort, obs_id)'; % SFirst values
        meanST(ids_sort, obs_id)'];
errs = [sdSF(ids_sort, obs_id)'; 
        sdST(ids_sort, obs_id)'];

b = bar(dat', ...
        'BarWidth', bw);
b(1).FaceColor = c1;
b(2).FaceColor = c2;

hold on

% Get x coordinates of bars
[ngroups, nbars] = size(dat);
x = nan(ngroups, nbars);
for i = 1:ngroups
    x(i,:) = b(i).XEndPoints;
end
% Plot the error bars
for i = 1:ngroups
    errorbar(x(i,:), dat(i,:), errs(i,:), 'k', 'linestyle', 'none');
end

grid on

legend(labs, 'fontsize', fleg)
%xlabel('parameter')
ylabel('Sobol indices')
set(gca, 'FontSize', fsize)
ylim([0,ymax])

xticks(1:length(ids_sort))
xticklabels(names(ids_sort))
temp = sprintf('Observable: Bone mineral density');
title(temp)


% AT1R
nexttile(2);
obs_id = 1;
% sort values by ST
[~, ids_st] = sort(meanST(ids, obs_id), 'descend');

ids_sort = ids(ids_st);

dat = [meanSF(ids_sort, obs_id)'; % SFirst values
        meanST(ids_sort, obs_id)'];
errs = [sdSF(ids_sort, obs_id)'; 
        sdST(ids_sort, obs_id)'];


b = bar(dat', ...
        'BarWidth', bw);
b(1).FaceColor = c1;
b(2).FaceColor = c2;

hold on

% Get x coordinates of bars
[ngroups, nbars] = size(dat);
x = nan(ngroups, nbars);
for i = 1:ngroups
    x(i,:) = b(i).XEndPoints;
end
% Plot the error bars
for i = 1:ngroups
    errorbar(x(i,:), dat(i,:), errs(i,:), 'k', 'linestyle', 'none');
end

grid on

legend(labs, 'fontsize', fleg)
ylabel('Sobol indices')
set(gca, 'FontSize', fsize)
ylim([0,ymax])

xticks(1:length(ids_sort))

xticklabels(names(ids_sort))
temp = sprintf('Observable: [AT1R-bound Ang II]');
title(temp)

f.Position = [100,100,width,height];

AddLetters2Plots(f, {'(A)', '(B)'},...
                'HShift', -0.06, 'VShift', -0.06,...
                'fontsize', 20)

%% Calciotropic hormones
%% additional figure
% % Determine which parameters have ST > threshold for 
% observations included
th = 0.02; 
obs_ids = [4,9];
ids = find(any(maxSF(:,obs_ids) > th, 2));

ymax = 0.4;
f=figure(5);
clf;
width = 1200;
height = 900;
f.Position = [100,100,width,height];

tiledlayout(2,1);


% PTH
nexttile(1);
obs_id = 4;
% sort values by ST
[~, ids_st] = sort(meanST(ids, obs_id), 'descend');

ids_sort = ids(ids_st);

dat = [meanSF(ids_sort, obs_id)'; % SFirst values
        meanST(ids_sort, obs_id)'];
errs = [sdSF(ids_sort, obs_id)'; 
        sdST(ids_sort, obs_id)'];

b = bar(dat', ...
        'BarWidth', bw);
b(1).FaceColor = c1;
b(2).FaceColor = c2;

hold on

% Get x coordinates of bars
[ngroups, nbars] = size(dat);
x = nan(ngroups, nbars);
for i = 1:ngroups
    x(i,:) = b(i).XEndPoints;
end
% Plot the error bars
for i = 1:ngroups
    errorbar(x(i,:), dat(i,:), errs(i,:), 'k', 'linestyle', 'none');
end

grid on
legend(labs, 'fontsize', fleg)
%xlabel('parameter')
ylabel('Sobol indices')
set(gca, 'FontSize', fsize)
ylim([0,ymax])

xticks(1:length(ids_sort))

xticklabels(names(ids_sort))
title('Observable: [PTH]')

legend(labs, 'fontsize', fleg)
%xlabel('parameter')
ylabel('Sobol indices')
set(gca, 'FontSize', fsize)
ylim([0,ymax])

xticks(1:length(ids_sort))

xticklabels(names(ids_sort))


% Ctriol
nexttile(2);
obs_id = 9;
% sort values by ST
[~, ids_st] = sort(meanST(ids, obs_id), 'descend');

ids_sort = ids(ids_st);

dat = [meanSF(ids_sort, obs_id)'; % SFirst values
        meanST(ids_sort, obs_id)'];
errs = [sdSF(ids_sort, obs_id)'; 
        sdST(ids_sort, obs_id)'];

b = bar(dat', ...
        'BarWidth', bw);
b(1).FaceColor = c1;
b(2).FaceColor = c2;

hold on


% Get x coordinates of bars
[ngroups, nbars] = size(dat);
x = nan(ngroups, nbars);
for i = 1:ngroups
    x(i,:) = b(i).XEndPoints;
end
% Plot the error bars
for i = 1:ngroups
    errorbar(x(i,:), dat(i,:), errs(i,:), 'k', 'linestyle', 'none');
end

grid on

legend(labs, 'fontsize', fleg)
%xlabel('parameter')
ylabel('Sobol indices')
%set(gcf, 'Position', [100, 100, 1000, 800]);
set(gca, 'FontSize', fsize)
ylim([0,ymax])

xticks(1:length(ids_sort))

xticklabels(names(ids_sort))
title('Observable: [Calcitriol]')

legend(labs, 'fontsize', fleg)
%xlabel('parameter')
ylabel('Sobol indices')
set(gca, 'FontSize', fsize)
ylim([0,ymax])

xticks(1:length(ids_sort))

xticklabels(names(ids_sort))

AddLetters2Plots(f, {'(A)', '(B)'},...
                'HShift', -0.06, 'VShift', -0.06,...
                'fontsize', 20)

%% additional figure
% % Determine which parameters have ST > threshold for 
% observations included
th = 0.05;
obs_ids = [3,4,6,9];
ids = find(any(maxSF(:,obs_ids) > th, 2));

ymax = 0.5;
f=figure(4);
clf;
width = 1200;
height = 900;
f.Position = [100,100,width,height];

tiledlayout(2,1);

% CaConc
nexttile(1);
obs_id = 3;
% sort values by ST
[~, ids_st] = sort(meanST(ids, obs_id), 'descend');

ids_sort = ids(ids_st);

dat = [meanSF(ids_sort, obs_id)'; % SFirst values
        meanST(ids_sort, obs_id)'];
errs = [sdSF(ids_sort, obs_id)'; 
        sdST(ids_sort, obs_id)'];

b = bar(dat', ...
        'BarWidth', bw);
b(1).FaceColor = c1;
b(2).FaceColor = c2;

hold on

% Get x coordinates of bars
[ngroups, nbars] = size(dat);
x = nan(ngroups, nbars);
for i = 1:ngroups
    x(i,:) = b(i).XEndPoints;
end
% Plot the error bars
for i = 1:ngroups
    errorbar(x(i,:), dat(i,:), errs(i,:), 'k', 'linestyle', 'none');
end

grid on

legend(labs, 'fontsize', fleg)
ylabel('Sobol indices')
set(gca, 'FontSize', fsize)
ylim([0,ymax])

xticks(1:length(ids_sort))

xticklabels(names(ids_sort))
temp = sprintf('Observable: %s', sRes.Observables{obs_id});
title('Observable: [Calcium]')
%title(temp)


% PRC
nexttile(2);
obs_id = 6;
% sort values by ST
[~, ids_st] = sort(meanST(ids, obs_id), 'descend');

ids_sort = ids(ids_st);

dat = [meanSF(ids_sort, obs_id)'; % SFirst values
        meanST(ids_sort, obs_id)'];
errs = [sdSF(ids_sort, obs_id)'; 
        sdST(ids_sort, obs_id)'];

b = bar(dat', ...
        'BarWidth', bw);
b(1).FaceColor = c1;
b(2).FaceColor = c2;

hold on

% Get x coordinates of bars
[ngroups, nbars] = size(dat);
x = nan(ngroups, nbars);
for i = 1:ngroups
    x(i,:) = b(i).XEndPoints;
end
% Plot the error bars
for i = 1:ngroups
    errorbar(x(i,:), dat(i,:), errs(i,:), 'k', 'linestyle', 'none');
end

grid on

legend(labs, 'fontsize', fleg)
%xlabel('parameter')
ylabel('Sobol indices')
set(gca, 'FontSize', fsize)
ylim([0,ymax])

xticks(1:length(ids_sort))

xticklabels(names(ids_sort))
title('Observable: [Renin]')




AddLetters2Plots(f, {'(A)', '(B)'},...
                'HShift', -0.06, 'VShift', -0.06,...
                'fontsize', 20)
