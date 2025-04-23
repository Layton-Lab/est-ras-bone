% Run ARB and ACEi at different ages

clear all;

tic
tvals = 20:0.25:80;

pct_ARB = 0.9359; % average
pct_ACEi = 0.956; % average 

% sim1: no RASi + EST def
fprintf('no RASi \n')
y1 = do_RASi_sim(0, 0, 0, 0, 0, 0, 1, tvals);

age_vals = [20, 30, 40, 50, 55, 60, 65, 70]; % ages to simulate
Y_ageARB = cell(size(age_vals));
Y_ageACEi = cell(size(age_vals));
for ii = 1:length(age_vals)
    age = age_vals(ii);
    disp(age);
    % ARB sim
    y = do_RASi_sim(1, 0, pct_ARB, 0, age, 0,1,tvals);
    Y_ageARB{ii} = y;
    % ACE i sim
    y = do_RASi_sim(0, 1, 0, pct_ACEi, 0, age, 1, tvals);
    Y_ageACEi{ii} = y;
end

%% Plot results
% fig specs
lw = 7;
cmap = parula(length(age_vals) + 2);
c1 = cmap(1,:);

xminmax=[20,80];
xlab = 't (years)';
fsize = 16;
labsARB = cell(length(age_vals) + 1, 1);
labsACEi = cell(length(age_vals) + 1, 1);

fleg = 12;

% plots
f = figure(1);
clf;
width = (2/3)*1600;
height = 900;
f.Position = [100,100,width,height];
tiledlayout(2,2);

% BMD (ARB)
nexttile(1);
id = 34; 
hold on
plot(tvals, y1(:,id)*100,'linewidth',lw,...
                    'color',c1)
labsARB{1} = 'no RASi';
for ii = 1:length(age_vals)
    labsARB{ii+1} = sprintf('%i years (ARB)', age_vals(ii));
    y = Y_ageARB{ii};
    c = cmap(ii + 1, :);
    if mod(ii, 3) == 1
        ls = ':';
    elseif mod(ii,3) == 2
        ls = '--';
    else
        ls = '-.';
    end
    plot(tvals,y(:,id)*100,'linewidth',lw,'color',c,'linestyle',ls)
end
grid on
xlim(xminmax)
ylim([50,120])
xlabel(xlab)
ylabel('Bone mineral density (%)')
set(gca, 'fontsize', fsize)
legend(labsARB,'location', 'southwest', 'fontsize', fleg)

% AT1R (ARB)
nexttile(2);
id = 32;
hold on
plot(tvals, y1(:,id),'linewidth',lw,...
                    'color',c1)
for ii = 1:length(age_vals)
    y = Y_ageARB{ii};
    c = cmap(ii + 1, :);
    if mod(ii, 3) == 1
        ls = ':';
    elseif mod(ii,3) == 2
        ls = '--';
    else
        ls = '-.';
    end
    plot(tvals,y(:,id),'linewidth',lw,'color',c,'linestyle',ls)
end
grid on
xlim(xminmax)
xlabel(xlab)
ylim([0,14])
ylabel('[AT1R-bound Ang II] (pmol/L)')
set(gca, 'fontsize', fsize)
legend(labsARB,'location', 'northwest', 'fontsize', fleg)

% BMD (ACEi)
nexttile(3);
id = 34; 
hold on
plot(tvals, y1(:,id)*100,'linewidth',lw,...
                    'color',c1)
labsACEi{1} = 'no RASi';
for ii = 1:length(age_vals)
    labsACEi{ii+1} = sprintf('%i years (ACEi)', age_vals(ii));
    y = Y_ageACEi{ii};
    c = cmap(ii + 1, :);
    if mod(ii, 3) == 1
        ls = ':';
    elseif mod(ii,3) == 2
        ls = '--';
    else
        ls = '-.';
    end
    plot(tvals,y(:,id)*100,'linewidth',lw,'color',c,'linestyle',ls)
end
grid on
xlim(xminmax)
ylim([50,120])
xlabel(xlab)
ylabel('Bone mineral density (%)')
set(gca, 'fontsize', fsize)
legend(labsACEi,'location', 'southwest', 'fontsize', fleg)

% AT1R (ACEi)
nexttile(4);
id = 32;
hold on
plot(tvals, y1(:,id),'linewidth',lw,...
                    'color',c1)
for ii = 1:length(age_vals)
    y = Y_ageACEi{ii};
    c = cmap(ii + 1, :);
    if mod(ii, 3) == 1
        ls = ':';
    elseif mod(ii,3) == 2
        ls = '--';
    else
        ls = '-.';
    end
    plot(tvals,y(:,id),'linewidth',lw,'color',c,'linestyle',ls)
end
grid on
xlim(xminmax)
xlabel(xlab)
ylim([0,14])
ylabel('[AT1R-bound Ang II] (pmol/L)')
set(gca, 'fontsize', fsize)
legend(labsACEi,'location', 'northwest', 'fontsize', fleg)

AddLetters2Plots(figure(1), {'(A1)', '(B1)', '(A2)', '(B2)'},...
                'HShift', -0.06, 'VShift', -0.06,...
                'fontsize', 20)





