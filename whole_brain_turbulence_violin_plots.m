% whole_brain_turbulence_violin_plots.m
%
% PURPOSE
% -------
% Planned whole-brain turbulence statistics and figures. This script loads
% turbulence_results.mat, extracts whole-brain turbulence at selected lambda
% values, tests predefined condition-group contrasts, applies FDR correction
% across lambdas within each contrast, plots paired violin figures, and saves
% the analysis table/arrays.
%
% INPUT FILE
% ----------
% turbulence_results.mat
%   Created by whole_brain_turbulence.m. Must contain output.Turbulence and
%   output.LAMBDA.
%
% REQUIRED HELPER FILES
% ---------------------
% Violin.m
%
% OUTPUT FILE
% -----------
% whole_brain_turbulence_lambda_results.mat
%   Saves target_lambdas, selected_lambdas, lambda_indices,
%   wholeBrainByLambda, and allResults.
%
% MAIN PARAMETERS TO CHANGE
% -------------------------
% target_lambdas:
%   Lambda values requested for statistics/plots.
%
% alpha:
%   FDR threshold.
%
% Condition indices:
%   MEMORY, COUNTING, J1, J4, J5, J7, J8.
%
% analyses:
%   Planned paired contrasts and direction labels.
%
% RUN ORDER
% ---------
% Run after:
%   whole_brain_turbulence.m
%
% NOTES
% -----
% This script uses whole-brain turbulence values only. It does not perform
% RSN/node-level analysis.

clear;
clc;
close all;

%% Settings
target_lambdas = [0.01 0.12 0.30];
alpha = 0.05;

MEMORY = 1;
COUNTING = 2;
J1 = 3;
J4 = 6;
J5 = 7;
J7 = 9;
J8 = 10;

analyses = { ...
    'Control vs Meditation', ...
    [MEMORY COUNTING], ...
    J1:J8, ...
    'Meditation - Control'; ...
    'Counting vs DeepMeditation', ...
    COUNTING, ...
    [J7 J8], ...
    'DeepMeditation - Counting'; ...
    'J1:J4 vs J5:J8', ...
    J1:J4, ...
    J5:J8, ...
    'J5:J8 - J1:J4'};

%% Load results
resultsFile = 'turbulence_results.mat';
loadedResults = load(resultsFile, 'output');

if ~isfield(loadedResults, 'output')
    error(['The file %s does not contain the variable output. ' ...
        'Rerun whole_brain_turbulence.m after its save call has been ' ...
        'updated to use MAT-file version 7.3.'], resultsFile);
end

output = loadedResults.output;

if ~isfield(output, 'Turbulence') || ~isfield(output, 'LAMBDA')
    error(['The variable output in %s is incomplete. Expected fields ' ...
        'output.Turbulence and output.LAMBDA. Rerun whole_brain_turbulence.m.'], ...
        resultsFile);
end

Turbulence = output.Turbulence;
LAMBDA = output.LAMBDA(:);

nSub = size(Turbulence, 1);
nCond = size(Turbulence, 2);
nLambda = length(target_lambdas);
nAnalyses = size(analyses, 1);

if nCond < 10
    error('Expected at least 10 conditions: Memory, Counting, J1-J8.');
end

%% Extract whole-brain turbulence for all requested lambdas
selected_lambdas = NaN(nLambda,1);
lambda_indices = NaN(nLambda,1);
wholeBrainByLambda = NaN(nSub, nCond, nLambda);

fprintf('\n========================================\n');
fprintf('WHOLE-BRAIN TURBULENCE ANALYSIS\n');
fprintf('========================================\n');
fprintf('Subjects in data = %d\n', nSub);
fprintf('Available lambda range = %.6f to %.6f\n', min(LAMBDA), max(LAMBDA));

fprintf('\nLambda selection\n');
fprintf('----------------------------------------\n');

for iLam = 1:nLambda
    target_lambda = target_lambdas(iLam);
    [~, lambda_index] = min(abs(LAMBDA - target_lambda));

    selected_lambdas(iLam) = LAMBDA(lambda_index);
    lambda_indices(iLam) = lambda_index;

    fprintf('Requested %.6f -> selected %.6f (index %d)\n', ...
        target_lambda, selected_lambdas(iLam), lambda_index);

    for sub = 1:nSub
        for cond = 1:nCond
            thisTurbulence = Turbulence{sub, cond};

            if isempty(thisTurbulence)
                continue
            end

            wholeBrainByLambda(sub, cond, iLam) = ...
                thisTurbulence(lambda_index);
        end
    end
end

%% Run each planned analysis
allResults = struct();

for iAnalysis = 1:nAnalyses
    analysisName = analyses{iAnalysis,1};
    condA = analyses{iAnalysis,2};
    condB = analyses{iAnalysis,3};
    diffLabel = analyses{iAnalysis,4};

    groupA = NaN(nSub, nLambda);
    groupB = NaN(nSub, nLambda);
    diffValues = NaN(nSub, nLambda);

    p_raw = NaN(nLambda,1);
    p_fdr = NaN(nLambda,1);
    effect_size = NaN(nLambda,1);
    n_pairs = NaN(nLambda,1);

    for iLam = 1:nLambda
        X = wholeBrainByLambda(:,:,iLam);

        groupA(:,iLam) = mean(X(:, condA), 2, 'omitnan');
        groupB(:,iLam) = mean(X(:, condB), 2, 'omitnan');
        diffValues(:,iLam) = groupB(:,iLam) - groupA(:,iLam);

        valid = ~isnan(groupA(:,iLam)) & ~isnan(groupB(:,iLam));
        x = groupA(valid,iLam);
        y = groupB(valid,iLam);

        n_pairs(iLam) = numel(x);

        if n_pairs(iLam) >= 2
            p_raw(iLam) = signrank(x, y);
            effect_size(iLam) = paired_rank_biserial(x, y);
        end
    end

    p_fdr = fdr_bh(p_raw, alpha);

    fprintf('\n========================================\n');
    fprintf('%s\n', analysisName);
    fprintf('Paired test direction: %s\n', diffLabel);
    fprintf('FDR correction across lambdas for this analysis\n');
    fprintf('Positive r means group B > group A\n');
    fprintf('========================================\n');

    for iLam = 1:nLambda
        fprintf('\nRequested lambda = %.6f | selected lambda = %.6f\n', ...
            target_lambdas(iLam), selected_lambdas(iLam));
        fprintf('Mean A = %.6f\n', mean(groupA(:,iLam), 'omitnan'));
        fprintf('Mean B = %.6f\n', mean(groupB(:,iLam), 'omitnan'));
        fprintf('Mean B-A = %.6f\n', mean(diffValues(:,iLam), 'omitnan'));
        fprintf('Pairs = %d\n', n_pairs(iLam));
        fprintf('Raw p = %.6g\n', p_raw(iLam));
        fprintf('FDR p = %.6g\n', p_fdr(iLam));
        fprintf('Rank-biserial r = %.4f\n', effect_size(iLam));

        fprintf('Subject values: A | B | B-A\n');
        validSubjects = find(~isnan(groupA(:,iLam)) & ~isnan(groupB(:,iLam)));
        for idx = 1:length(validSubjects)
            sub = validSubjects(idx);
            fprintf('  sub %02d: %.6f | %.6f | %.6f\n', ...
                sub, groupA(sub,iLam), groupB(sub,iLam), diffValues(sub,iLam));
        end
    end

    allResults(iAnalysis).name = analysisName;
    allResults(iAnalysis).condA = condA;
    allResults(iAnalysis).condB = condB;
    allResults(iAnalysis).target_lambdas = target_lambdas;
    allResults(iAnalysis).selected_lambdas = selected_lambdas;
    allResults(iAnalysis).lambda_indices = lambda_indices;
    allResults(iAnalysis).groupA = groupA;
    allResults(iAnalysis).groupB = groupB;
    allResults(iAnalysis).diffValues = diffValues;
    allResults(iAnalysis).p_raw = p_raw;
    allResults(iAnalysis).p_fdr = p_fdr;
    allResults(iAnalysis).effect_size = effect_size;
    allResults(iAnalysis).n_pairs = n_pairs;

    plot_lambda_difference_violins( ...
        groupA, ...
        groupB, ...
        target_lambdas, ...
        selected_lambdas, ...
        p_raw, ...
        p_fdr, ...
        analysisName, ...
        analyses{iAnalysis,1});
end

%% Save analysis outputs
save('whole_brain_turbulence_lambda_results.mat', ...
    'target_lambdas', ...
    'selected_lambdas', ...
    'lambda_indices', ...
    'wholeBrainByLambda', ...
    'allResults');

fprintf('\nSaved: whole_brain_turbulence_lambda_results.mat\n');

%% Local helper functions
function plot_lambda_difference_violins(groupA, groupB, target_lambdas, ...
    selected_lambdas, p_raw, p_fdr, analysisName, comparisonLabel)

    figure('Color','w','Position',[100 100 1100 650]);
    hold on;

    colorA = [0.20 0.35 0.70];
    colorB = [0.75 0.25 0.25];

    nLambda = size(groupA,2);
    xCenters = 1:nLambda;
    xA = xCenters - 0.18;
    xB = xCenters + 0.18;

    for iLam = 1:nLambda
        valuesA = groupA(:,iLam);
        valuesB = groupB(:,iLam);
        valuesA = valuesA(~isnan(valuesA));
        valuesB = valuesB(~isnan(valuesB));

        Violin({valuesA}, xA(iLam), ...
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

        Violin({valuesB}, xB(iLam), ...
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

    for sub = 1:size(groupA,1)
        for iLam = 1:nLambda
            a = groupA(sub,iLam);
            b = groupB(sub,iLam);

            if ~isnan(a) && ~isnan(b)
                plot([xA(iLam) xB(iLam)], [a b], '-', ...
                    'Color', [0.78 0.78 0.78], ...
                    'LineWidth', 0.8);
            end

            if ~isnan(a)
                scatter(xA(iLam), a, 34, ...
                    'MarkerFaceColor', colorA, ...
                    'MarkerEdgeColor', 'w', ...
                    'LineWidth', 0.7, ...
                    'MarkerFaceAlpha', 0.85);
            end

            if ~isnan(b)
                scatter(xB(iLam), b, 34, ...
                    'MarkerFaceColor', colorB, ...
                    'MarkerEdgeColor', 'w', ...
                    'LineWidth', 0.7, ...
                    'MarkerFaceAlpha', 0.85);
            end
        end
    end

    for iLam = 1:nLambda
        medianA = median(groupA(:,iLam), 'omitnan');
        medianB = median(groupB(:,iLam), 'omitnan');

        plot([xA(iLam) - 0.09, xA(iLam) + 0.09], ...
            [medianA, medianA], ...
            'k-', 'LineWidth', 2.4);

        plot([xB(iLam) - 0.09, xB(iLam) + 0.09], ...
            [medianB, medianB], ...
            'k-', 'LineWidth', 2.4);
    end

    allValues = [groupA(:); groupB(:)];
    yMin = min(allValues, [], 'omitnan');
    yMax = max(allValues, [], 'omitnan');
    yRange = yMax - yMin;

    if yRange == 0
        yRange = max(abs(yMax), 1) * 0.1;
    end

    dataPad = 0.08 * yRange;
    labelY = yMax + 0.04 * yRange;
    stepY = 0.06 * yRange;
    labelTop = labelY + (nLambda-1)*stepY + 0.08*yRange;

    for iLam = 1:nLambda
        text(xCenters(iLam), labelY + (iLam-1)*stepY, ...
            sprintf('raw p=%.3g\nFDR p=%.3g', p_raw(iLam), p_fdr(iLam)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontName','Arial', ...
            'FontSize',10);
    end

    xLabels = cell(nLambda,1);
    for iLam = 1:nLambda
        xLabels{iLam} = sprintf('\\lambda = %.2f', target_lambdas(iLam));
    end

    xticks(1:nLambda);
    xticklabels(xLabels);
    ylabel('Whole-brain turbulence');
    legendLabels = strsplit(comparisonLabel, ' vs ');
    if numel(legendLabels) ~= 2
        legendLabels = {'Group A','Group B'};
    end

    hA = plot(NaN, NaN, 'o', ...
        'MarkerFaceColor', colorA, ...
        'MarkerEdgeColor', colorA, ...
        'LineStyle', 'none');
    hB = plot(NaN, NaN, 'o', ...
        'MarkerFaceColor', colorB, ...
        'MarkerEdgeColor', colorB, ...
        'LineStyle', 'none');

    title(sprintf('%s | Whole-brain turbulence by lambda', analysisName));
    legend([hA hB], legendLabels, 'Location','bestoutside');

    xlim([0.45 nLambda+0.55]);
    ylim([yMin - dataPad, labelTop]);

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

function p_adj = fdr_bh(p, alpha)
% Benjamini-Hochberg FDR adjusted p-values.

    if nargin < 2
        alpha = 0.05;
    end

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

    %#ok<NASGU>
end
