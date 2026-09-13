% hopf_capacity_susceptibility_violin_plots.m
%
% PURPOSE
% Plot subject-level Hopf information capacity and susceptibility results.
% The script creates violin plots for all conditions from
% hopf_subjectwise_IC_SUS.mat.
%
% INPUT FILE
% hopf_subjectwise_IC_SUS.mat
%   Created by hopf_capacity_susceptibility.m. Must contain InfoCapacity,
%   Susceptibility, condition_order, and OptimalG.
%
% REQUIRED HELPER FILES
% Violin.m
%
% MAIN PARAMETERS TO CHANGE
% cols:
%   Plot color map.
%
% RUN ORDER
% Run after:
%   hopf_capacity_susceptibility.m

clear; close all; clc;

results_file = 'hopf_subjectwise_IC_SUS.mat';
if ~isfile(results_file)
    error('Cannot find %s. Run hopf_capacity_susceptibility.m first.', ...
        results_file);
end

load(results_file);
% must contain:
% InfoCapacity (nSub x nCond)
% Susceptibility (nSub x nCond)
% condition_order
% OptimalG

IC  = InfoCapacity;
SUS = Susceptibility;

nCond = size(IC,2);
nSub  = size(IC,1);

cols = lines(nCond);
% FIGURE 1 — ALL CONDITIONS (IC)

figure('Color','w','Position',[100 100 1300 550]);
hold on;

for c = 1:nCond

    data = IC(:,c);
    data = data(~isnan(data));

    % VIOLIN
    Violin({data}, c, ...
        'ViolinColor',{cols(c,:)}, ...
        'ViolinAlpha',{0.28}, ...
        'EdgeColor',[0.2 0.2 0.2], ...
        'ShowData',false, ...
        'ShowBox',false, ...
        'ShowWhiskers',false, ...
        'ShowMedian',false, ...
        'Width',0.32);

    % SUBJECT DOTS (each dot = subject mean across simulations)
    jitter = 0.05;
    scatter(c + jitter*randn(numel(data),1), data, ...
        35, ...
        'filled', ...
        'MarkerFaceColor',cols(c,:), ...
        'MarkerEdgeColor','w');

    % MEDIAN
    med = median(data,'omitnan');
    plot([c-0.18 c+0.18], [med med], 'k-', 'LineWidth',2);

    % MEAN
    mu = mean(data,'omitnan');
    scatter(c, mu, 80, 'd', 'filled', ...
        'MarkerFaceColor',[0.85 0 0], ...
        'MarkerEdgeColor','k');
end

set(gca, ...
    'XTick',1:nCond, ...
    'XTickLabel',condition_order, ...
    'FontSize',11, ...
    'LineWidth',1.2, ...
    'Box','off');

ylabel('Information Capacity');
title('All Conditions — Subject-Level Distributions');
% FIGURE 2 — ALL CONDITIONS (SUS)

figure('Color','w','Position',[100 100 1300 550]);
hold on;

for c = 1:nCond

    data = SUS(:,c);
    data = data(~isnan(data));

    Violin({data}, c, ...
        'ViolinColor',{cols(c,:)}, ...
        'ViolinAlpha',{0.28}, ...
        'EdgeColor',[0.2 0.2 0.2], ...
        'ShowData',false, ...
        'ShowBox',false, ...
        'ShowWhiskers',false, ...
        'ShowMedian',false, ...
        'Width',0.32);

    jitter = 0.05;
    scatter(c + jitter*randn(numel(data),1), data, ...
        35, 'filled', ...
        'MarkerFaceColor',cols(c,:), ...
        'MarkerEdgeColor','w');

    med = median(data,'omitnan');
    plot([c-0.18 c+0.18], [med med], 'k-', 'LineWidth',2);

    mu = mean(data,'omitnan');
    scatter(c, mu, 80, 'd', 'filled', ...
        'MarkerFaceColor',[0.85 0 0], ...
        'MarkerEdgeColor','k');
end

set(gca, ...
    'XTick',1:nCond, ...
    'XTickLabel',condition_order, ...
    'FontSize',11, ...
    'LineWidth',1.2, ...
    'Box','off');

ylabel('Susceptibility');
title('All Conditions — Subject-Level Distributions');
% PLANNED COMPARISONS (NO FDR, PAIRED TESTS)

MEMORY = 1;
COUNTING = 2;
J1 = 3;
J4 = 6;
J5 = 7;
J7 = 9;
J8 = 10;

comparisons = {
    'Control vs Meditation',        [MEMORY COUNTING], J1:J8;
    'Counting vs DeepMeditation',   COUNTING,          [J7 J8];
    'J1:J4 vs J5:J8',              J1:J4,             J5:J8;
};
% FUNCTION FOR COMPARISON PLOTS

plot_comparison = @(A,B,labels,titleStr,ylabelStr) ...
    comparison_violin(A,B,labels,titleStr,ylabelStr,cols);
% RUN COMPARISONS — IC

for i = 1:size(comparisons,1)

    name = comparisons{i,1};

    Aidx = comparisons{i,2};

    Bidx = comparisons{i,3};

    A = mean(IC(:,Aidx),2,'omitnan');
    B = mean(IC(:,Bidx),2,'omitnan');

    valid = ~isnan(A) & ~isnan(B);
    A = A(valid);
    B = B(valid);

    p = signrank(A,B);

    comparison_violin(A,B,{'A','B'},name,'IC',cols);

    text(1.5, max([A;B])*1.05, sprintf('p = %.3g',p), ...
        'HorizontalAlignment','center', ...
        'FontSize',12);
end
% RUN COMPARISONS — SUS

for i = 1:size(comparisons,1)

    name = comparisons{i,1};

    Aidx = comparisons{i,2};

    Bidx = comparisons{i,3};

    A = mean(SUS(:,Aidx),2,'omitnan');
    B = mean(SUS(:,Bidx),2,'omitnan');

    valid = ~isnan(A) & ~isnan(B);
    A = A(valid);
    B = B(valid);

    p = signrank(A,B);

    comparison_violin(A,B,{'A','B'},name,'SUS',cols);

    text(1.5, max([A;B])*1.05, sprintf('p = %.3g',p), ...
        'HorizontalAlignment','center', ...
        'FontSize',12);
end
% SUPPORT FUNCTION

function comparison_violin(A,B,labels,titleStr,ylabelStr,cols)

figure('Color','w','Position',[200 200 500 450]);
hold on;

data = {A,B};

for i = 1:2

    Violin({data{i}}, i, ...
        'ViolinColor',{cols(i,:)}, ...
        'ViolinAlpha',{0.25}, ...
        'EdgeColor',[0.2 0.2 0.2], ...
        'ShowData',false, ...
        'Width',0.35);

    jitter = 0.06;
    scatter(i + jitter*randn(numel(data{i}),1), data{i}, ...
        40, 'filled', ...
        'MarkerFaceColor',cols(i,:), ...
        'MarkerEdgeColor','w');

    med = median(data{i},'omitnan');
    plot([i-0.15 i+0.15],[med med],'k-','LineWidth',2);
end

set(gca, ...
    'XTick',1:2, ...
    'XTickLabel',labels, ...
    'FontSize',11, ...
    'LineWidth',1.2, ...
    'Box','off');

ylabel(ylabelStr);
title(titleStr);
end
