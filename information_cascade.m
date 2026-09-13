% information_cascade.m
%
% PURPOSE
% Compute information cascade from whole-brain turbulence results. The script
% loads local Kuramoto synchronization values from turbulence_results.mat,
% estimates transfer between adjacent lambda scales over time, summarizes
% planned condition contrasts, plots violin comparisons, and saves the
% subject-level results/statistics.
%
% METHOD SUMMARY
% output.R{sub,cond} is interpreted as lambda x parcels x time. For each
% subject and condition, the script correlates local synchronization at
% lambda(k+1), time t+1 with lambda(k), time t across parcels. Significant
% correlations with p < 0.05 are kept, absolute values are averaged per
% lambda transition, then averaged across lambda transitions.
%
% INPUT FILE
% turbulence_results.mat
%   Created by whole_brain_turbulence.m. Must contain output.R and
%   output.LAMBDA.
%
% REQUIRED HELPER FILES
% Violin.m
%
% OUTPUT FILE
% information_cascade_results.mat
%   Saves LAMBDA, TransferLambda_sub, InformationCascade_sub, group summary
%   vectors, comparison definitions, p-values, effect sizes, and pair counts.
%
% MAIN PARAMETERS TO CHANGE
% conditionLabels:
%   Display names for the 10 conditions.
%
% Condition indices:
%   MEMORY, COUNTING, J1, ..., J8.
%
% comparisons:
%   Planned paired Wilcoxon signed-rank contrasts.
%
% Significance threshold:
%   The cascade calculation currently keeps correlations with pp < 0.05.
%
% RUN ORDER
% Run after:
%   whole_brain_turbulence.m
%
% Optional follow-up:
%   information_cascade_violin_plots.m can recreate the violin figures from
%   information_cascade_results.mat without recomputing the cascade.

clear;
clc;

%% Load turbulence output
input_file = 'turbulence_results.mat';
if ~isfile(input_file)
    error('Cannot find %s. Run whole_brain_turbulence.m first.', input_file);
end

load(input_file, 'output');

LAMBDA = output.LAMBDA(:);
R_all = output.R;

NLAMBDA = length(LAMBDA);
NSUB = size(R_all,1);
NCOND = size(R_all,2);

conditionLabels = { ...
    'Memory', 'Counting', ...
    'J1', 'J2', 'J3', 'J4', 'J5', 'J6', 'J7', 'J8'};

%% Storage
% TransferLambda_sub(k,sub,cond) corresponds to transfer into lambda k.
% The first lambda has no previous scale, so it remains NaN.
TransferLambda_sub = NaN(NLAMBDA, NSUB, NCOND);
InformationCascade_sub = NaN(NSUB, NCOND);

%% Main computation
for sub = 1:NSUB
    fprintf('Subject %d / %d\n', sub, NSUB);

    for cond = 1:NCOND
        R = R_all{sub,cond};

        if isempty(R)
            continue
        end

        R = fix_R_dimensions(R, NLAMBDA);

        for ilam = 1:NLAMBDA-1
            % Original method:
            % corr(enstrophy(lambda+1,:,2:end)', enstrophy(lambda,:,1:end-1)')
            nextScale = squeeze(R(ilam+1,:,2:end))';
            thisScale = squeeze(R(ilam,:,1:end-1))';

            [cc, pp] = corr(nextScale, thisScale, 'Rows','pairwise');
            sig = pp(:) < 0.05;

            if any(sig)
                TransferLambda_sub(ilam+1,sub,cond) = ...
                    mean(abs(cc(sig)), 'omitnan');
            end
        end

        InformationCascade_sub(sub,cond) = ...
            mean(TransferLambda_sub(2:NLAMBDA,sub,cond), 'omitnan');
    end
end

%% Group summaries used in the rest of this project
MEMORY = 1;
COUNTING = 2;
J1 = 3;
J2 = 4;
J3 = 5;
J4 = 6;
J5 = 7;
J6 = 8;
J7 = 9;
J8 = 10;

Control = mean(InformationCascade_sub(:, [MEMORY COUNTING]), 2, 'omitnan');
Counting = InformationCascade_sub(:, COUNTING);
Meditation = mean(InformationCascade_sub(:, J1:J8), 2, 'omitnan');
DeepMeditation = mean(InformationCascade_sub(:, [J7 J8]), 2, 'omitnan');
EarlyMeditation = mean(InformationCascade_sub(:,J1:J4),2,'omitnan');
LateMeditation  = mean(InformationCascade_sub(:,J5:J8),2,'omitnan');

fprintf('\n========================================\n');
fprintf('INFORMATION CASCADE RESULTS\n');
fprintf('========================================\n');
fprintf('Lambda range: %.6f to %.6f\n', min(LAMBDA), max(LAMBDA));
fprintf('Number of lambdas: %d\n', NLAMBDA);
fprintf('Subjects: %d\n', NSUB);
fprintf('Conditions: %d\n', NCOND);

fprintf('\nCondition means across subjects\n');
fprintf('----------------------------------------\n');
for cond = 1:min(NCOND, length(conditionLabels))
    fprintf('%s = %.6f\n', ...
        conditionLabels{cond}, ...
        mean(InformationCascade_sub(:,cond), 'omitnan'));
end

fprintf('\nProject group means across subjects\n');
fprintf('----------------------------------------\n');
fprintf('Control          = %.6f\n', mean(Control, 'omitnan'));
fprintf('Counting         = %.6f\n', mean(Counting, 'omitnan'));
fprintf('Meditation       = %.6f\n', mean(Meditation, 'omitnan'));
fprintf('DeepMeditation   = %.6f\n', mean(DeepMeditation, 'omitnan'));

%% Paired statistics
comparisons = { ...
    'Control vs Meditation',      Control,         Meditation; ...
    'Counting vs DeepMeditation', Counting,        DeepMeditation; ...
    'J1-J4 vs J5-J8',             EarlyMeditation, LateMeditation};

nTests = size(comparisons,1);
p_raw = NaN(nTests,1);
effect_size = NaN(nTests,1);
n_pairs = NaN(nTests,1);

for iTest = 1:nTests
    x = comparisons{iTest,2};
    y = comparisons{iTest,3};
    valid = ~isnan(x) & ~isnan(y);

    x = x(valid);
    y = y(valid);
    n_pairs(iTest) = numel(x);

    if n_pairs(iTest) >= 2
        p_raw(iTest) = signrank(x, y);
        effect_size(iTest) = paired_rank_biserial(x, y);
    end
end



fprintf('\nPaired Wilcoxon tests on information cascade\n');
fprintf('Positive r means second condition > first condition\n');
fprintf('----------------------------------------\n');
for iTest = 1:nTests
    fprintf('%s\n', comparisons{iTest,1});
    fprintf('  pairs        = %d\n', n_pairs(iTest));
    fprintf('  raw p        = %.6g\n', p_raw(iTest));
    fprintf('  rank-biserial r = %.4f\n', effect_size(iTest));
end

% FIGURE 1 : CONTROL VS MEDITATION

figure('Color','w','Position',[100 100 650 600]);
hold on;

plotData = [Control Meditation];
plotLabels = {'Control','Meditation'};

colors = [
    0.20 0.35 0.70
    0.20 0.55 0.35
];

plot_violin_comparison(plotData,plotLabels,colors,...
    sprintf('Control vs Meditation\np = %.3g',p_raw(1)));

% FIGURE 2 : COUNTING VS DEEP MEDITATION

figure('Color','w','Position',[100 100 650 600]);
hold on;

plotData = [Counting DeepMeditation];
plotLabels = {'Counting','Deep Meditation'};

colors = [
    0.30 0.45 0.80
    0.75 0.25 0.25
];

plot_violin_comparison(plotData,plotLabels,colors,...
    sprintf('Counting vs Deep Meditation\np = %.3g',p_raw(2)));

% FIGURE 3 : J1-J4 VS J5-J8

figure('Color','w','Position',[100 100 650 600]);
hold on;

plotData = [EarlyMeditation LateMeditation];
plotLabels = {'J1-J4','J5-J8'};

colors = [
    0.35 0.60 0.85
    0.85 0.45 0.25
];

plot_violin_comparison(plotData,plotLabels,colors,...
    sprintf('J1-J4 vs J5-J8\np = %.3g',p_raw(3)));

%% Save
save('information_cascade_results.mat', ...
    'LAMBDA', ...
    'TransferLambda_sub', ...
    'InformationCascade_sub', ...
    'Control', ...
    'Counting', ...
    'Meditation', ...
    'DeepMeditation', ...
    'EarlyMeditation', ...
    'LateMeditation', ...
    'conditionLabels', ...
    'comparisons', ...
    'p_raw', ...
    'effect_size', ...
    'n_pairs');

fprintf('\nSaved: information_cascade_results.mat\n');

%% Local helper
function R = fix_R_dimensions(R, nLambda)
% Return R as lambda x parcels x time.

    dims = size(R);
    lambda_dim = find(dims == nLambda, 1);
    parcel_dim = find(dims == 1000, 1);
    if isempty(parcel_dim)
        parcel_dim = find(dims == 200, 1);
    end

    if isempty(lambda_dim)
        error('Cannot find lambda dimension in output.R cell.');
    end

    if isempty(parcel_dim)
        error('Cannot find parcel dimension in output.R cell.');
    end

    time_dim = setdiff(1:ndims(R), [lambda_dim parcel_dim]);

    if numel(time_dim) ~= 1
        error('Expected one time dimension in output.R cell.');
    end

    R = permute(R, [lambda_dim parcel_dim time_dim]);
end


function r = paired_rank_biserial(x, y)
% Positive values mean y is larger than x.

    d = y(:) - x(:);
    d = d(d ~= 0 & ~isnan(d));

    if isempty(d)
        r = 0;
        return
    end

    ranks = tiedrank(abs(d));
    wPos = sum(ranks(d > 0));
    wNeg = sum(ranks(d < 0));
    r = (wPos - wNeg) / (wPos + wNeg);
end

function add_pvalue_bar(x1, x2, y, labelText)
% Draw comparison bracket and p-value label.

    yl = ylim;
    tickHeight = 0.015 * (max(yl) - min(yl));
    plot([x1 x1 x2 x2], [y y+tickHeight y+tickHeight y], ...
        'k-', 'LineWidth', 1.2);
    text(mean([x1 x2]), y + 1.5*tickHeight, labelText, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 11, ...
        'FontName', 'Arial');
end

function plot_violin_comparison(plotData,plotLabels,colors,plotTitle)

NSUB = size(plotData,1);

for iGroup = 1:2

    values = plotData(:,iGroup);
    values = values(~isnan(values));

    Violin({values}, iGroup, ...
        'ViolinColor',{colors(iGroup,:)}, ...
        'ViolinAlpha',{0.28}, ...
        'EdgeColor',[0.25 0.25 0.25], ...
        'BoxColor',[0.15 0.15 0.15], ...
        'MedianColor',[1 1 1], ...
        'ShowData',false, ...
        'ShowBox',false, ...
        'ShowWhiskers',false, ...
        'ShowMedian',true, ...
        'Width',0.32);

end

xJitter = [-0.03 0.03];

for sub = 1:NSUB

    y = plotData(sub,:);
    valid = ~isnan(y);

    if sum(valid)==2
        plot([1 2],y,'-', ...
            'Color',[0.75 0.75 0.75], ...
            'LineWidth',0.8);
    end

    for iGroup = 1:2

        if ~isnan(y(iGroup))

            scatter(iGroup+xJitter(iGroup), ...
                y(iGroup), ...
                38, ...
                'MarkerFaceColor',colors(iGroup,:), ...
                'MarkerEdgeColor','w', ...
                'LineWidth',0.7, ...
                'MarkerFaceAlpha',0.85);

        end
    end
end

medians = median(plotData,1,'omitnan');

for iGroup = 1:2

    plot([iGroup-0.18 iGroup+0.18], ...
         [medians(iGroup) medians(iGroup)], ...
         'k-', ...
         'LineWidth',2.4);

end

yMin = min(plotData(:),[],'omitnan');
yMax = max(plotData(:),[],'omitnan');
yRange = yMax-yMin;

if yRange==0
    yRange = 1;
end

barY = yMax + 0.08*yRange;

add_pvalue_bar(1,2,barY,'');

xlim([0.45 2.55]);
ylim([yMin-0.08*yRange yMax+0.20*yRange]);

xticks([1 2]);
xticklabels(plotLabels);

ylabel('Information Cascade');
title(plotTitle);

set(gca,...
    'FontName','Arial',...
    'FontSize',13,...
    'LineWidth',1.2,...
    'TickDir','out',...
    'Box','off');

grid on;
ax = gca;
ax.GridAlpha = 0.12;
ax.XGrid = 'off';

end
