
% hopf_model_fit_plots.m
%
% PURPOSE
% Create diagnostic figures from Hopf SC fit results: optimal-G violins,
% Control vs Meditation optimal-G comparison, representative empirical vs
% simulated FC matrices, and group-level fit error as a function of G.
%
% INPUT FILE
% fit_results_file
%   Hopf result file created by hopf_model_fit.m or by a merged/extended G
%   search. Must contain Cfg, G_range, OptimalG, FitErrorMean, EmpFC, Fdiff,
%   Tmax, SC, subjects_clean, and condition_order.
%
% REQUIRED HELPER FILES
% Violin.m
%
% OUTPUT FILES
% Figure1_OptimalG_Violins.png
% Figure1B_OptimalG_Control_vs_Meditation.png
% Figure2_Empirical_vs_Simulated_FC.png
% Figure3_GroupFit_Error_vs_G.png
%
% MAIN PARAMETERS TO CHANGE
% fit_results_file:
%   Which Hopf result file to plot.
%
% condition_order:
%   Condition labels shown in plots. The script uses labels from the fit file
%   if available; otherwise it falls back to the listed defaults.
%
% control_idx / meditation_idx:
%   Condition groups for the planned optimal-G comparison.
%
% cond:
%   Representative condition used for empirical-vs-simulated FC.
%
% RUN ORDER
% Run after:
%   hopf_model_fit.m

clear;
clc;
close all;

fit_results_file = 'hopf_sc_fit_results_schaefer1000_G0_to_4p5.mat';
if ~isfile(fit_results_file)
    error('Cannot find %s. Run hopf_model_fit.m first or change fit_results_file.', ...
        fit_results_file);
end
load(fit_results_file);

% CONDITION NAMES
if ~exist('condition_order', 'var')
    condition_order = { ...
        'Memory', ...
        'Counting', ...
        'J1', ...
        'J2', ...
        'J3', ...
        'J4', ...
        'J5', ...
        'J6', ...
        'J7', ...
        'J8'};
end

% FILTER

fnq = 1/(2*Cfg.TR);

Wn = [Cfg.filt.lb/fnq ...
      Cfg.filt.ub/fnq];

[bfilt,afilt] = butter(2,Wn);

% FIGURE 1
% VIOLIN PLOTS OF OPTIMAL G

figure('Color','w','Position',[100 100 1200 600]);
hold on;

nCond = size(OptimalG,2);
cols = lines(nCond);

for c = 1:nCond

    data = OptimalG(:,c);
    data = data(~isnan(data));

    % --- VIOLIN 
    Violin({data}, c, ...
        'ViolinColor',{cols(c,:)}, ...
        'ViolinAlpha',{0.28}, ...
        'EdgeColor',[0.25 0.25 0.25], ...
        'ShowData',false, ...
        'ShowBox',false, ...
        'ShowWhiskers',false, ...
        'ShowMedian',false, ...
        'Width',0.32);

    % --- SUBJECT DOTS (jittered, structured)
    jitter = 0.04;

    scatter(c + jitter*randn(numel(data),1), ...
            data, ...
            38, ...
            'MarkerFaceColor',cols(c,:), ...
            'MarkerEdgeColor','w', ...
            'LineWidth',0.7, ...
            'MarkerFaceAlpha',0.85);

    % --- MEDIAN
    med = median(data);

    plot([c-0.18 c+0.18], [med med], ...
        'k-', ...
        'LineWidth',2.4);

    % --- MEAN
    mu = mean(data,'omitnan');

    scatter(c, mu, ...
        90, ...
        'd', ...
        'filled', ...
        'MarkerFaceColor',[0.85 0 0], ...
        'MarkerEdgeColor','k', ...
        'LineWidth',1.0);
    
end

set(gca,...
    'XTick',1:nCond,...
    'XTickLabel',condition_order,...
    'FontSize',12,...
    'LineWidth',1.2,...
    'TickDir','out',...
    'Box','off');

ylabel('Optimal G');
xlabel('Condition');
title('Distribution of Optimal G Across Subjects');

grid on;
ax = gca;
ax.GridAlpha = 0.12;
ax.XGrid = 'off';

saveas(gcf,'Figure1_OptimalG_Violins.png');

% FIGURE 1B
% CONTROL vs MEDITATION OPTIMAL G

control_idx = [1 2];     % Memory + Counting
meditation_idx = 3:10;   % J1-J8

control_G = mean(OptimalG(:,control_idx),2,'omitnan');
meditation_G = mean(OptimalG(:,meditation_idx),2,'omitnan');

valid = ~isnan(control_G) & ~isnan(meditation_G);

control_clean = control_G(valid);
meditation_clean = meditation_G(valid);

if numel(control_clean) > 1
    p_control_meditation = signrank(control_clean,meditation_clean);
else
    p_control_meditation = NaN;
end

fprintf('\n');
fprintf('Control vs Meditation Optimal G:\n');
fprintf('Subjects used = %d\n',numel(control_clean));
fprintf('Control mean = %.5f\n',mean(control_clean,'omitnan'));
fprintf('Meditation mean = %.5f\n',mean(meditation_clean,'omitnan'));
fprintf('Wilcoxon signed-rank p = %.6f\n',p_control_meditation);

figure('Color','w','Position',[150 150 650 520]);
hold on;

comparison_data = {control_clean,meditation_clean};
comparison_cols = [0.25 0.45 0.85; 0.85 0.25 0.25];

for iGroup = 1:2

    data = comparison_data{iGroup};

    Violin({data}, iGroup, ...
        'ViolinColor',{comparison_cols(iGroup,:)}, ...
        'ViolinAlpha',{0.28}, ...
        'EdgeColor',[0.25 0.25 0.25], ...
        'ShowData',false, ...
        'ShowBox',false, ...
        'ShowWhiskers',false, ...
        'ShowMedian',false, ...
        'Width',0.35);

    jitter = 0.05;

    scatter(iGroup + jitter*randn(numel(data),1), ...
            data, ...
            42, ...
            'MarkerFaceColor',comparison_cols(iGroup,:), ...
            'MarkerEdgeColor','w', ...
            'LineWidth',0.7, ...
            'MarkerFaceAlpha',0.85);

    med = median(data,'omitnan');

    plot([iGroup-0.18 iGroup+0.18], [med med], ...
        'k-', ...
        'LineWidth',2.4);

    mu = mean(data,'omitnan');

    scatter(iGroup, mu, ...
        90, ...
        'd', ...
        'filled', ...
        'MarkerFaceColor',[0.85 0 0], ...
        'MarkerEdgeColor','k', ...
        'LineWidth',1.0);

end

for sub = 1:numel(control_clean)
    plot([1 2], ...
         [control_clean(sub) meditation_clean(sub)], ...
         '-', ...
         'Color',[0.75 0.75 0.75], ...
         'LineWidth',0.8);
end

set(gca,...
    'XTick',1:2,...
    'XTickLabel',{'Control','Meditation'},...
    'FontSize',12,...
    'LineWidth',1.2,...
    'TickDir','out',...
    'Box','off');

ylabel('Optimal G');
xlabel('Condition group');
title('Optimal G: Control vs Meditation');

grid on;
ax = gca;
ax.GridAlpha = 0.12;
ax.XGrid = 'off';

yVals = [control_clean; meditation_clean];
yVals = yVals(~isnan(yVals));

if isempty(yVals)
    yMin = 0;
    yMax = 1;
else
    yMin = min(yVals);
    yMax = max(yVals);
end

yRange = yMax - yMin;
if yRange == 0
    yRange = max(abs(yMax),1) * 0.1;
end

bracketY = yMax + 0.08*yRange;
textY = yMax + 0.14*yRange;
ylim([yMin - 0.08*yRange, yMax + 0.24*yRange]);

plot([1 1 2 2], ...
     [bracketY bracketY+0.025*yRange bracketY+0.025*yRange bracketY], ...
     'k-', ...
     'LineWidth',1.2);

text(1.5, textY, sprintf('p = %.3g',p_control_meditation), ...
    'HorizontalAlignment','center', ...
    'FontSize',12);

if p_control_meditation < 0.05
    if p_control_meditation < 0.001
        sig_label = '***';
    elseif p_control_meditation < 0.01
        sig_label = '**';
    else
        sig_label = '*';
    end

    text(1.5, textY + 0.07*yRange, sig_label, ...
        'HorizontalAlignment','center', ...
        'FontSize',18, ...
        'FontWeight','bold');
end

saveas(gcf,'Figure1B_OptimalG_Control_vs_Meditation.png');

% FIGURE 2 EMPIRICAL VS SIMULATED FC


cond = 1;

groupMeanG = mean(OptimalG(:,cond),'omitnan');

[~,sub] = min(abs( ...
    OptimalG(:,cond) - groupMeanG));

fprintf('\n');
fprintf('Representative subject selected:\n');
fprintf('Subject = %s\n',string(subjects_clean{sub}));
fprintf('Condition = %s\n',string(condition_order{cond}));

G = OptimalG(sub,cond);

fprintf('Optimal G = %.3f\n',G);

f_diff  = Fdiff{sub,cond};
thisTmax = Tmax(sub,cond);

%% IMPORTANT:
% Same stochastic setup used during fitting

rng(1);

sim_ts = simulate_hopf_sc( ...
    SC,...
    G,...
    f_diff,...
    thisTmax,...
    Cfg);

sim_filt = filter_bold( ...
    sim_ts,...
    bfilt,...
    afilt,...
    0);

simFC = corrcoef(sim_filt');

%% Plot FC matrices

figure( ...
    'Color','w', ...
    'Position',[100 100 1000 450]);

subplot(1,2,1)

imagesc(EmpFC{sub,cond});

axis square
colorbar

caxis([-1 1])

title(sprintf( ...
    'Empirical FC\n%s | %s', ...
    string(subjects_clean{sub}), ...
    string(condition_order{cond})));

subplot(1,2,2)

imagesc(simFC);

axis square
colorbar

caxis([-1 1])

title(sprintf( ...
    'Simulated FC\nG = %.2f', ...
    G));

colormap(parula)

saveas(gcf,'Figure2_Empirical_vs_Simulated_FC.png');
% FIGURE 3
% FIT ERROR VS G (MEAN ± SUBJECT VARIABILITY)
% FIGURE 3
% GROUP-LEVEL FIT: ERROR VS G

figure( ...
    'Color','w', ...
    'Position',[100 100 1200 700]);

hold on

cols = lines(nCond);

for cond = 1:nCond

    errmat = squeeze(FitErrorMean(:,cond,:));

    meanErr = mean(errmat,1,'omitnan');
    sdErr   = std(errmat,0,1,'omitnan');

    fill( ...
        [G_range fliplr(G_range)], ...
        [meanErr+sdErr fliplr(meanErr-sdErr)], ...
        cols(cond,:), ...
        'FaceAlpha',0.12, ...
        'EdgeColor','none', ...
        'HandleVisibility','off');

    plot( ...
        G_range, ...
        meanErr, ...
        '-o', ...
        'Color',cols(cond,:), ...
        'LineWidth',2.2, ...
        'MarkerSize',7, ...
        'MarkerFaceColor','none');

end

xlabel('G (global coupling)','FontSize',16);
ylabel('Fit error','FontSize',16);

title('Group-level fit: error vs G', ...
    'FontSize',18, ...
    'FontWeight','bold');

legend(condition_order, ...
    'Location','eastoutside');

set(gca,...
    'FontSize',14,...
    'LineWidth',1.2,...
    'Box','on');

grid on

saveas(gcf,'Figure3_GroupFit_Error_vs_G.png');

fprintf('\n');
fprintf('All figures saved successfully.\n');
% LOCAL FUNCTIONS

function BOLD_filt = filter_bold(BOLD,bfilt,afilt,excTp)

nAreas = size(BOLD,1);
nTime = size(BOLD,2);

BOLD_filt = zeros(nAreas,nTime-2*excTp);

for seed = 1:nAreas

    sig = BOLD(seed,:) - mean(BOLD(seed,:));

    sig = detrend(sig);

    sig = filtfilt(bfilt,afilt,sig);

    if excTp > 0
        sig = sig(excTp+1:end-excTp);
    end

    BOLD_filt(seed,:) = sig;

end

end

function xs = simulate_hopf_sc(C,G,f_diff,Tmax,Cfg)

nAreas = Cfg.nNodes;

TR = Cfg.TR;

dt = 0.1*TR/2;

omega = repmat(2*pi*f_diff(:),1,2);

omega(:,1) = -omega(:,1);

dsig = sqrt(dt)*Cfg.sig;

wC = G*C;

sumC = repmat(sum(wC,2),1,2);

a = Cfg.a*ones(nAreas,2);

z = 0.1*ones(nAreas,2);

for t = 0:dt:2000

    z = hopf_step( ...
        z,a,omega,wC,sumC,dt,dsig);

end

xs = zeros(Tmax,nAreas);

nn = 0;

for t = 0:dt:((Tmax-1)*TR)

    z = hopf_step( ...
        z,a,omega,wC,sumC,dt,dsig);

    if abs(mod(t,TR)) < 0.01

        nn = nn + 1;

        xs(nn,:) = z(:,1)';

    end

end

xs = xs';

end

function z = hopf_step(z,a,omega,wC,sumC,dt,dsig)

suma = wC*z - sumC.*z;

zz = z(:,end:-1:1);

z = z + dt * ( ...
    a.*z + ...
    zz.*omega - ...
    z.*(z.*z + zz.*zz) + ...
    suma ) + ...
    dsig * randn(size(z));

end
