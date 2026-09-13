% metastability_kuramoto.m
%
% PURPOSE
% -------
% Compute global Kuramoto order-parameter summaries from Hilbert phases, then
% run subject-level metastability statistics and plots.
%
% The first part loads phases created by BOLD2hilbertupdated.m, computes:
%   OP_all : Kuramoto order-parameter time series per subject-condition
%   Sync   : mean order parameter per subject-condition
%   Meta   : standard deviation of order parameter per subject-condition
% and appends them to preprocessed_data.mat.
%
% The second part compares metastability across conditions using Friedman and
% paired Wilcoxon signed-rank tests, creates violin plots, and saves summary
% statistics.
%
% INPUT FILES
% -----------
% preprocessed_data.mat
%   Must contain Phase, created by BOLD2hilbertupdated.m.
%
% REQUIRED HELPER FILES
% ---------------------
% metastability_kuramoto_helper.m
% Violin.m
%
% OUTPUT FILES
% ------------
% preprocessed_data.mat
%   Appends OP_all, Sync, and Meta.
%
% metastability_statistics.mat
%   Saves Meta, p-values, FDR-corrected pairwise p-values, comparison pairs,
%   and rank-biserial effect sizes for planned contrasts.
%
% MAIN PARAMETERS TO CHANGE
% -------------------------
% Condition indices:
%   MEMORY, COUNTING, J1, ..., J8
%   These assume condition order:
%   Memory, Counting, J1, J2, J3, J4, J5, J6, J7, J8.
%
% Planned contrasts:
%   control          = mean(Memory, Counting)
%   meditation       = mean(J1:J8)
%   early meditation = mean(J1:J4)
%   late meditation  = mean(J5:J8)
%   deep meditation  = mean(J7,J8)
%
% RUN ORDER
% ---------
% Run after:
%   BOLD2hilbertupdated.m
%
% Run before, if needed:
%   downstream scripts that use preprocessed_data.mat with Meta/Sync/OP_all.

clear; clc;

input_file = 'preprocessed_data.mat';
if ~isfile(input_file)
    error('Cannot find %s. Run BOLD2hilbertupdated.m first.', input_file);
end

load(input_file)

nAreas = size(Phase{1,1}, 1);

OP_all = cell(size(Phase));
Sync = nan(size(Phase));  
Meta = nan(size(Phase));

for i = 1:size(Phase,1)
    for j = 1:size(Phase,2)

        if isempty(Phase{i,j})
            continue
        end

        Phase_BOLD = Phase{i,j};

        [OP, Synchro, Metasta] = metastability_kuramoto_helper(Phase_BOLD, nAreas);

        OP_all{i,j} = OP; %for subject per session you have the full time series so subject x session where each cell is a full time series (vector)  
        Sync(i,j) = Synchro; %subject per session matrix
        Meta(i,j) = Metasta; %subj per session matrix 

    end
end

save('preprocessed_data.mat', 'OP_all', 'Meta', 'Sync', '-append')


% SUBJECT-LEVEL METASTABILITY ANALYSIS

fprintf('\n========================================\n');
fprintf('SUBJECT-LEVEL METASTABILITY ANALYSIS\n');
fprintf('========================================\n');

%% Condition indices

MEMORY   = 1;
COUNTING = 2;

J1 = 3;
J2 = 4;
J3 = 5;
J4 = 6;

J5 = 7;
J6 = 8;
J7 = 9;
J8 = 10;

% FIGURE 1: ALL CONDITIONS

figure('Color','w','Position',[100 100 1200 600]);
hold on

for cond = 1:10

    values = Meta(:,cond);
    values = values(~isnan(values));

    Violin({values}, cond, ...
        'ViolinAlpha',{0.30}, ...
        'ShowData',true, ...
        'ShowMedian',true, ...
        'ShowBox',false, ...
        'ShowWhiskers',false);

end

xticks(1:10)
xticklabels({'Memory','Counting','J1','J2','J3','J4','J5','J6','J7','J8'})

ylabel('Metastability')
title('Metastability across all conditions')

[p_friedman,tbl,stats] = friedman(Meta,1,'off');

fprintf('\nALL CONDITIONS\n');
fprintf('Friedman p = %.6g\n',p_friedman);

%% Pairwise post-hoc tests

pairs = nchoosek(1:10,2);

pairwise_p = nan(size(pairs,1),1);

for k = 1:size(pairs,1)

    c1 = pairs(k,1);
    c2 = pairs(k,2);

    pairwise_p(k) = signrank( ...
        Meta(:,c1), ...
        Meta(:,c2));

end

pairwise_p_fdr = fdr_bh(pairwise_p);

% FIGURE 2: CONTROL VS MEDITATION

control = mean(Meta(:,[MEMORY COUNTING]),2);

meditation = mean(Meta(:,J1:J8),2);

p_control_meditation = signrank(control,meditation);

r_control_meditation = ...
    paired_rank_biserial(control,meditation);

figure('Color','w','Position',[100 100 700 600]);
hold on

Violin({control},1,...
    'ViolinColor',{[0.2 0.4 0.8]},...
    'ViolinAlpha',{0.3},...
    'ShowData',true);

Violin({meditation},2,...
    'ViolinColor',{[0.8 0.3 0.3]},...
    'ViolinAlpha',{0.3},...
    'ShowData',true);

for s = 1:length(control)

    if ~isnan(control(s)) && ~isnan(meditation(s))

        plot([1 2], ...
             [control(s) meditation(s)], ...
             '-', ...
             'Color',[0.8 0.8 0.8]);

    end

end

xticks([1 2])
xticklabels({'Control','Meditation'})

ylabel('Metastability')

title(sprintf( ...
    'Control vs Meditation\np = %.4f | r = %.3f', ...
    p_control_meditation, ...
    r_control_meditation));

% FIGURE 3: J1:J4 VS J5:J8

early = mean(Meta(:,J1:J4),2);

late = mean(Meta(:,J5:J8),2);

p_early_late = signrank(early,late);

r_early_late = paired_rank_biserial(early,late);

figure('Color','w','Position',[100 100 700 600]);
hold on

Violin({early},1,...
    'ViolinColor',{[0.2 0.4 0.8]},...
    'ViolinAlpha',{0.3},...
    'ShowData',true);

Violin({late},2,...
    'ViolinColor',{[0.8 0.3 0.3]},...
    'ViolinAlpha',{0.3},...
    'ShowData',true);

for s = 1:length(early)

    if ~isnan(early(s)) && ~isnan(late(s))

        plot([1 2], ...
             [early(s) late(s)], ...
             '-', ...
             'Color',[0.8 0.8 0.8]);

    end

end

xticks([1 2])
xticklabels({'J1-J4','J5-J8'})

ylabel('Metastability')

title(sprintf( ...
    'J1:J4 vs J5:J8\np = %.4f | r = %.3f', ...
    p_early_late, ...
    r_early_late));

% FIGURE 4: COUNTING VS DEEP MEDITATION

counting = Meta(:,COUNTING);

deepMeditation = mean(Meta(:,[J7 J8]),2);

p_counting_deep = signrank(counting,deepMeditation);

r_counting_deep = paired_rank_biserial( ...
    counting, ...
    deepMeditation);

figure('Color','w','Position',[100 100 700 600]);
hold on

Violin({counting},1,...
    'ViolinColor',{[0.2 0.4 0.8]},...
    'ViolinAlpha',{0.3},...
    'ShowData',true);

Violin({deepMeditation},2,...
    'ViolinColor',{[0.8 0.3 0.3]},...
    'ViolinAlpha',{0.3},...
    'ShowData',true);

for s = 1:length(counting)

    if ~isnan(counting(s)) && ~isnan(deepMeditation(s))

        plot([1 2], ...
             [counting(s) deepMeditation(s)], ...
             '-', ...
             'Color',[0.8 0.8 0.8]);

    end

end

xticks([1 2])
xticklabels({'Counting','Deep Meditation'})

ylabel('Metastability')

title(sprintf( ...
    'Counting vs Deep Meditation\np = %.4f | r = %.3f', ...
    p_counting_deep, ...
    r_counting_deep));

% SAVE RESULTS

save('metastability_statistics.mat', ...
    'Meta', ...
    'p_friedman', ...
    'pairwise_p', ...
    'pairwise_p_fdr', ...
    'pairs', ...
    'p_control_meditation', ...
    'p_early_late', ...
    'p_counting_deep', ...
    'r_control_meditation', ...
    'r_early_late', ...
    'r_counting_deep');

fprintf('\nSaved: metastability_statistics.mat\n');

% HELPER FUNCTIONS

function r = paired_rank_biserial(x,y)

d = y(:) - x(:);

d = d(~isnan(d));
d = d(d~=0);

if isempty(d)
    r = NaN;
    return
end

ranks = tiedrank(abs(d));

wPos = sum(ranks(d>0));
wNeg = sum(ranks(d<0));

r = (wPos - wNeg) / (wPos + wNeg);

end

function p_adj = fdr_bh(p)

p_adj = nan(size(p));

valid = ~isnan(p);

p_valid = p(valid);

[p_sorted,idx] = sort(p_valid);

m = numel(p_sorted);

adj = p_sorted .* m ./ (1:m)';

adj = min(adj,1);

for k = m-1:-1:1
    adj(k) = min(adj(k),adj(k+1));
end

tmp = nan(m,1);
tmp(idx) = adj;

p_adj(valid) = tmp;

end
