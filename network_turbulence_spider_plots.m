% network_turbulence_spider_plots.m
%
% PURPOSE
% Network-level turbulence plots and statistics. This script loads
% turbulence_results.mat, extracts local synchronization output.R at a
% selected lambda, converts it into node-level turbulence, averages parcels
% within Schaefer-7 resting-state networks, and produces spider/violin plots
% for planned condition-group comparisons.
%
% INPUT FILES
% turbulence_results.mat
%   Created by whole_brain_turbulence.m. Must contain output.R and
%   output.LAMBDA.
%
% schaefer1000_7N_coords_decolab.csv
%   Used for Schaefer network labels.
%
% REQUIRED HELPER FILES
% load_schaefer_coords.m
% Violin.m
%
% OUTPUT
% Figures and printed statistics are produced in MATLAB.
%
% MAIN PARAMETERS TO CHANGE
% target_lambda:
%   Lambda used for network-level plots/statistics. The closest available
%   lambda in output.LAMBDA is selected automatically.
%
% run_spider calls:
%   Planned condition contrasts and plot labels.
%
% RUN ORDER
% Run after:
%   whole_brain_turbulence.m

clear;
clc;
% LOAD DATA

input_file = 'turbulence_results.mat';
if ~isfile(input_file)
    error('Cannot find %s. Run whole_brain_turbulence.m first.', input_file);
end

load(input_file);

% Select the spatial scale used for the spider plots.
% The script automatically uses the closest available lambda in output.LAMBDA.
target_lambda = 0.12;
LAMBDA = output.LAMBDA(:);
[~, lambda_index] = min(abs(LAMBDA - target_lambda));
selected_lambda = LAMBDA(lambda_index);

fprintf('\n========================================\n');
fprintf('SPIDER PLOT LAMBDA SELECTION\n');
fprintf('Requested lambda = %.4f\n', target_lambda);
fprintf('Using lambda     = %.4f\n', selected_lambda);
fprintf('Lambda index     = %d\n', lambda_index);
fprintf('========================================\n');

[~, Labels] = load_schaefer_coords('schaefer1000_7N_coords_decolab.csv');
NPARCELS = length(Labels);
% DEFINE NETWORKS

Networks = cell(NPARCELS,1);

for i = 1:NPARCELS
    label = Labels{i};

    if contains(label,'Vis')
        Networks{i} = 'VIS';
    elseif contains(label,'SomMot')
        Networks{i} = 'SOM';
    elseif contains(label,'DorsAttn')
        Networks{i} = 'DAN';
    elseif contains(label,'SalVentAttn')
        Networks{i} = 'VEN';
    elseif contains(label,'Limbic')
        Networks{i} = 'LIM';
    elseif contains(label,'Cont')
        Networks{i} = 'FPN';
    elseif contains(label,'Default')
        Networks{i} = 'DMN';
    else
        Networks{i} = 'UNK';
    end
end

NetNames = {'VIS','SOM','DAN','VEN','LIM','FPN','DMN'};
Nnet = length(NetNames);
% RUN ANALYSES

run_spider(output, Networks, NetNames, Nnet, ...
    [1 2], ...
    [3 4 5 6 7 8 9 10], ...
    'Control vs Meditation', ...
    lambda_index, length(LAMBDA), selected_lambda);

run_spider(output, Networks, NetNames, Nnet, ...
    [2], ...
    [9 10], ...
    'Counting vs DeepMeditation', ...
    lambda_index, length(LAMBDA), selected_lambda);

run_spider(output, Networks, NetNames, Nnet, ...
    [3 4 5 6], ...
    [7 8 9 10], ...
    'J1:J4 vs J5:J8', ...
    lambda_index, length(LAMBDA), selected_lambda);
% MAIN FUNCTION

function run_spider(output, Networks, NetNames, Nnet, condA, condB, ...
    title_str, lambda_index, nLambda, selected_lambda)

NSUB = size(output.R,1);
NPARCELS = length(Networks);

NodeA = NaN(NPARCELS, NSUB);
NodeB = NaN(NPARCELS, NSUB);

for sub = 1:NSUB

    NodeA(:,sub) = compute_node_metastability( ...
        output, sub, condA, lambda_index, nLambda, NPARCELS);

    NodeB(:,sub) = compute_node_metastability( ...
        output, sub, condB, lambda_index, nLambda, NPARCELS);
end
% PAPER-STYLE NODE DIFFERENCE AND RSN COUNT

MeanA = mean(NodeA,2,'omitnan');
MeanB = mean(NodeB,2,'omitnan');

NodeDiff = abs(MeanB - MeanA);

% Select the top 15% largest node-level absolute differences.
threshold = quantile(NodeDiff, 0.85);
TopNodes = NodeDiff >= threshold;

Counts = zeros(Nnet,1);
RSN_A = NaN(1,Nnet);
RSN_B = NaN(1,Nnet);
p_raw = NaN(1,Nnet);
effect_size = NaN(1,Nnet);
n_pairs = NaN(1,Nnet);
RSN_SubA = NaN(NSUB,Nnet);
RSN_SubB = NaN(NSUB,Nnet);

for n = 1:Nnet
    idx = strcmp(Networks, NetNames{n});
    Counts(n) = sum(TopNodes(idx));
    RSN_A(n) = mean(MeanA(idx), 'omitnan');
    RSN_B(n) = mean(MeanB(idx), 'omitnan');

    subjA = mean(NodeA(idx,:), 1, 'omitnan')';
    subjB = mean(NodeB(idx,:), 1, 'omitnan')';
    RSN_SubA(:,n) = subjA;
    RSN_SubB(:,n) = subjB;

    valid = ~isnan(subjA) & ~isnan(subjB);

    n_pairs(n) = sum(valid);

    if n_pairs(n) >= 2
        p_raw(n) = signrank(subjA(valid), subjB(valid));
        effect_size(n) = paired_rank_biserial(subjA(valid), subjB(valid));
    end
end

p_fdr = fdr_bh(p_raw);

fprintf('\n========================================\n');
fprintf('%s | lambda = %.4f\n', title_str, selected_lambda);
fprintf('Node-level metastability difference, top 15%% counted by RSN\n');
fprintf('Top 15%% threshold = %.6f\n', threshold);
fprintf('----------------------------------------\n');
for n = 1:Nnet
    fprintf('%s | top nodes = %d\n', NetNames{n}, Counts(n));
end

fprintf('\nRSN turbulence statistics for paired violin plot\n');
fprintf('Paired Wilcoxon signrank tests across subjects, FDR across RSNs\n');
fprintf('Positive r means group B > group A\n');
fprintf('----------------------------------------\n');
for n = 1:Nnet
    fprintf(['%s | A = %.6f | B = %.6f | pairs = %d | ' ...
        'raw p = %.6g | FDR p = %.6g | r = %.4f\n'], ...
        NetNames{n}, RSN_A(n), RSN_B(n), n_pairs(n), ...
        p_raw(n), p_fdr(n), effect_size(n));
end
% PLOT

figure;

Counts_plot = [Counts; Counts(1)];
theta = linspace(0,2*pi,length(Counts_plot));

polarplot(theta, Counts_plot, 'LineWidth', 2);

ax = gca;
ax.ThetaTick = linspace(0,360,8);
ax.ThetaTickLabel = [NetNames NetNames(1)];

title(sprintf('%s | Top 15%% node differences by RSN | \\lambda = %.4f', ...
    title_str, selected_lambda));

legendLabels = strsplit(title_str, ' vs ');
if numel(legendLabels) ~= 2
    legendLabels = {'Group A','Group B'};
end

plot_rsn_turbulence_violins(RSN_SubA, RSN_SubB, NetNames, ...
    legendLabels, p_raw, p_fdr, title_str, selected_lambda);

if strcmp(title_str, 'Control vs Meditation')
    plot_control_meditation_spiders(RSN_A, RSN_B, NetNames, selected_lambda);
end

end
% RSN TURBULENCE VIOLIN PLOT

function plot_rsn_turbulence_violins(RSN_SubA, RSN_SubB, NetNames, ...
    legendLabels, p_raw, p_fdr, title_str, selected_lambda)

Nnet = length(NetNames);

figure('Color','w','Position',[100 100 1250 650]);
hold on;

colorA = [0.20 0.35 0.70];
colorB = [0.75 0.25 0.25];
xCenters = 1:Nnet;
xA = xCenters - 0.18;
xB = xCenters + 0.18;

for n = 1:Nnet
    a = RSN_SubA(:,n);
    b = RSN_SubB(:,n);
    a = a(~isnan(a));
    b = b(~isnan(b));

    Violin({a}, xA(n), ...
        'ViolinColor', {colorA}, ...
        'ViolinAlpha', {0.28}, ...
        'EdgeColor', [0.25 0.25 0.25], ...
        'BoxColor', [0.15 0.15 0.15], ...
        'MedianColor', [1 1 1], ...
        'ShowData', false, ...
        'ShowBox', false, ...
        'ShowWhiskers', false, ...
        'ShowMedian', true, ...
        'Width', 0.20);

    Violin({b}, xB(n), ...
        'ViolinColor', {colorB}, ...
        'ViolinAlpha', {0.28}, ...
        'EdgeColor', [0.25 0.25 0.25], ...
        'BoxColor', [0.15 0.15 0.15], ...
        'MedianColor', [1 1 1], ...
        'ShowData', false, ...
        'ShowBox', false, ...
        'ShowWhiskers', false, ...
        'ShowMedian', true, ...
        'Width', 0.20);
end

for sub = 1:size(RSN_SubA,1)
    for n = 1:Nnet
        a = RSN_SubA(sub,n);
        b = RSN_SubB(sub,n);

        if ~isnan(a) && ~isnan(b)
            plot([xA(n) xB(n)], [a b], '-', ...
                'Color', [0.78 0.78 0.78], ...
                'LineWidth', 0.8);
        end

        if ~isnan(a)
            scatter(xA(n), a, 30, ...
                'MarkerFaceColor', colorA, ...
                'MarkerEdgeColor', 'w', ...
                'LineWidth', 0.6, ...
                'MarkerFaceAlpha', 0.85);
        end

        if ~isnan(b)
            scatter(xB(n), b, 30, ...
                'MarkerFaceColor', colorB, ...
                'MarkerEdgeColor', 'w', ...
                'LineWidth', 0.6, ...
                'MarkerFaceAlpha', 0.85);
        end
    end
end

for n = 1:Nnet
    medianA = median(RSN_SubA(:,n), 'omitnan');
    medianB = median(RSN_SubB(:,n), 'omitnan');

    plot([xA(n)-0.08 xA(n)+0.08], [medianA medianA], ...
        'k-', 'LineWidth', 2.2);
    plot([xB(n)-0.08 xB(n)+0.08], [medianB medianB], ...
        'k-', 'LineWidth', 2.2);
end

allValues = [RSN_SubA(:); RSN_SubB(:)];
yMin = min(allValues, [], 'omitnan');
yMax = max(allValues, [], 'omitnan');
yRange = yMax - yMin;
if yRange == 0
    yRange = max(abs(yMax), 1) * 0.1;
end

labelY = yMax + 0.04*yRange;
labelStep = 0.035*yRange;

for n = 1:Nnet
    text(xCenters(n), labelY + mod(n-1,2)*labelStep, ...
        sprintf('p=%.3g\nFDR=%.3g', p_raw(n), p_fdr(n)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontName','Arial', ...
        'FontSize',9);
end

hA = plot(NaN, NaN, 'o', ...
    'MarkerFaceColor', colorA, ...
    'MarkerEdgeColor', colorA, ...
    'LineStyle','none');
hB = plot(NaN, NaN, 'o', ...
    'MarkerFaceColor', colorB, ...
    'MarkerEdgeColor', colorB, ...
    'LineStyle','none');

xlim([0.45 Nnet+0.55]);
ylim([yMin - 0.08*yRange, yMax + 0.20*yRange]);
xticks(1:Nnet);
xticklabels(NetNames);
ylabel('RSN node-level turbulence');
title(sprintf('%s | RSN turbulence at \\lambda = %.4f', ...
    title_str, selected_lambda));
legend([hA hB], legendLabels, 'Location','bestoutside');

set(gca, ...
    'FontName','Arial', ...
    'FontSize',13, ...
    'LineWidth',1.2, ...
    'TickDir','out', ...
    'Box','off');

grid on;
ax = gca;
ax.GridAlpha = 0.12;
ax.XGrid = 'off';

end
% EXTRA CONTROL/MEDITATION RSN SPIDER PLOTS

function plot_control_meditation_spiders(RSN_Control, RSN_Meditation, ...
    NetNames, selected_lambda)

plot_single_rsn_spider( ...
    RSN_Control, ...
    NetNames, ...
    [0.20 0.35 0.70], ...
    sprintf(['Control RSN turbulence | \\lambda = %.4f\n' ...
    'Control = mean(Memory, Counting)'], selected_lambda));

plot_single_rsn_spider( ...
    RSN_Meditation, ...
    NetNames, ...
    [0.75 0.25 0.25], ...
    sprintf(['Meditation RSN turbulence | \\lambda = %.4f\n' ...
    'Meditation = mean(J1:J8)'], selected_lambda));

figure('Color','w','Position',[100 100 850 700]);
theta = linspace(0, 2*pi, length(NetNames)+1);
polarplot(theta, [RSN_Control RSN_Control(1)], ...
    'LineWidth', 2.2, ...
    'Color', [0.20 0.35 0.70]);
hold on;
polarplot(theta, [RSN_Meditation RSN_Meditation(1)], ...
    'LineWidth', 2.2, ...
    'Color', [0.75 0.25 0.25]);
format_rsn_polar_axis(NetNames);
legend({'Control','Meditation'}, 'Location','bestoutside');
title(sprintf(['Control vs Meditation overlapping RSN turbulence | \\lambda = %.4f\n' ...
    'Lines show mean RSN node-level turbulence'], selected_lambda));

diffValues = RSN_Meditation - RSN_Control;
plot_single_rsn_spider( ...
    abs(diffValues), ...
    NetNames, ...
    [0.25 0.25 0.25], ...
    sprintf(['Meditation - Control RSN turbulence difference | \\lambda = %.4f\n' ...
    'Radius shows absolute magnitude; all signed values are Meditation - Control'], ...
    selected_lambda));

fprintf('\nControl/Meditation spider summary values\n');
fprintf('----------------------------------------\n');
for n = 1:length(NetNames)
    fprintf('%s | Control = %.6f | Meditation = %.6f | Meditation-Control = %.6f\n', ...
        NetNames{n}, RSN_Control(n), RSN_Meditation(n), diffValues(n));
end

end

function plot_single_rsn_spider(values, NetNames, color, titleText)

figure('Color','w','Position',[100 100 850 700]);
theta = linspace(0, 2*pi, length(values)+1);
polarplot(theta, [values values(1)], ...
    'LineWidth', 2.2, ...
    'Color', color);
format_rsn_polar_axis(NetNames);
title(titleText);

end

function format_rsn_polar_axis(NetNames)

ax = gca;
ax.ThetaTick = linspace(0,360,length(NetNames)+1);
ax.ThetaTickLabel = [NetNames NetNames(1)];
ax.FontName = 'Arial';
ax.FontSize = 13;
grid on;

end
% NODE-LEVEL METASTABILITY FUNCTION

function nodeValues = compute_node_metastability(output, sub, conditions, ...
    lambda_index, nLambda, nParcels)

nodeByCondition = NaN(nParcels, length(conditions));

for iCond = 1:length(conditions)

    c = conditions(iCond);
    X = output.R{sub,c};

    if isempty(X)
        continue
    end

    % After fixing, X has dimensions parcels x time for the selected lambda.
    X = fix_data(X, lambda_index, nLambda, nParcels);

    if size(X,1) ~= nParcels
        error('Unexpected parcel count after fixing output.R data.');
    end

    % Node-level metastability: standard deviation across time of the
    % local Kuramoto order parameter for each node.
    nodeByCondition(:,iCond) = std(X, 0, 2, 'omitnan');
end

% If a group contains multiple conditions, average node-level metastability
% across those conditions within subject.
nodeValues = mean(nodeByCondition, 2, 'omitnan');

end
% SAFE DATA FIX FUNCTION

function X = fix_data(X, lambda_index, nLambda, nParcels)

dims = size(X);
pdim = find(dims == nParcels, 1);

if isempty(pdim)
    error('Cannot find %d-parcel dimension', nParcels);
end

% If a lambda dimension exists, keep only the selected lambda.
lambda_dim = [];

if ndims(X) > 2
    lambda_dim = find(size(X) == nLambda & (1:ndims(X)) ~= pdim, 1);
end

if ~isempty(lambda_dim)
    idx = repmat({':'}, 1, ndims(X));
    idx{lambda_dim} = lambda_index;
    X = X(idx{:});

    % Recompute the parcel dimension after slicing.
    pdim = find(size(X) == nParcels, 1);
    if isempty(pdim)
        error('Cannot find %d-parcel dimension after selecting lambda.', nParcels);
    end
end

% Move parcels to rows.
X = permute(X, [pdim setdiff(1:ndims(X), pdim)]);

% Keep all remaining dimensions as observations.
% For selected-lambda output.R this gives parcels x time.
X = reshape(X, size(X,1), []);

end
% FDR CORRECTION

function p_adj = fdr_bh(p)

p_adj = NaN(size(p));
valid = ~isnan(p);
p_valid = p(valid);

if isempty(p_valid)
    return
end

[p_sorted, sortIdx] = sort(p_valid(:));
m = numel(p_sorted);

adj_sorted = p_sorted .* m ./ (1:m)';
adj_sorted = min(adj_sorted, 1);

for i = m-1:-1:1
    adj_sorted(i) = min(adj_sorted(i), adj_sorted(i+1));
end

tmp = NaN(m,1);
tmp(sortIdx) = adj_sorted;
p_adj(valid) = tmp;

end
% PAIRED RANK-BISERIAL EFFECT SIZE

function r = paired_rank_biserial(x, y)

% Positive values mean group B is larger than group A.
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
