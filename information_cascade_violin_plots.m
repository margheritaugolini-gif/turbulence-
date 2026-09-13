% information_cascade_violin_plots.m
%
% PURPOSE
% Recreate information-cascade violin plots from saved results without
% recomputing the cascade measure.
%
% INPUT FILE
% information_cascade_results.mat
%   Created by information_cascade.m.
%
% REQUIRED HELPER FILES
% Violin.m
%
% MAIN PARAMETERS TO CHANGE
% Colors and plotLabels in each figure block.
%
% RUN ORDER
% Run after:
%   information_cascade.m

clear;
clc;

results_file = 'information_cascade_results.mat';
if ~isfile(results_file)
    error('Cannot find %s. Run information_cascade.m first.', results_file);
end

load(results_file, ...
    'Control', ...
    'Counting', ...
    'Meditation', ...
    'DeepMeditation', ...
    'EarlyMeditation', ...
    'LateMeditation', ...
    'p_raw');

figure('Color','w','Position',[100 100 650 600]);
hold on;

plotData = [Control Meditation];
plotLabels = {'Control','Meditation'};

colors = [
    0.20 0.35 0.70
    0.20 0.55 0.35
];

plot_violin_comparison_local(plotData,plotLabels,colors,...
    sprintf('Control vs Meditation\np = %.3g',p_raw(1)));

figure('Color','w','Position',[100 100 650 600]);
hold on;

plotData = [Counting DeepMeditation];
plotLabels = {'Counting','Deep Meditation'};

colors = [
    0.30 0.45 0.80
    0.75 0.25 0.25
];

plot_violin_comparison_local(plotData,plotLabels,colors,...
    sprintf('Counting vs Deep Meditation\np = %.3g',p_raw(2)));

figure('Color','w','Position',[100 100 650 600]);
hold on;

plotData = [EarlyMeditation LateMeditation];
plotLabels = {'J1-J4','J5-J8'};

colors = [
    0.35 0.60 0.85
    0.85 0.45 0.25
];

plot_violin_comparison_local(plotData,plotLabels,colors,...
    sprintf('J1-J4 vs J5-J8\np = %.3g',p_raw(3)));

function add_pvalue_bar_local(x1, x2, y, labelText)
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

function plot_violin_comparison_local(plotData,plotLabels,colors,plotTitle)

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

add_pvalue_bar_local(1,2,barY,'');

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
