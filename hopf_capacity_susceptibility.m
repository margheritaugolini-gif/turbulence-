% hopf_capacity_susceptibility.m
%
% PURPOSE
% Compute subject-level information capacity and susceptibility from fitted
% SC-coupled Hopf models. For each subject and condition, the script uses the
% fitted OptimalG and empirical node frequencies, runs baseline and perturbed
% Hopf simulations, computes SC-weighted enstrophy, and summarizes the
% perturbation response.
%
% INPUT FILE
% fit_results_file
%   Hopf fit results created by hopf_model_fit.m or by a merged/extended G
%   search. Must contain OptimalG, Fdiff, Tmax, subjects_clean, and
%   condition_order. If SC is absent, the script loads it with
%   load_schaefer1000_sc.m.
%
% REQUIRED HELPER FILES
% load_schaefer1000_sc.m
% readNPY.m, only if loading SC from .npy
%
% OUTPUT FILE
% hopf_subjectwise_IC_SUS.mat
%   Saves InfoCapacity, Susceptibility, OptimalG, subjects_clean,
%   condition_order, and measure_weighting.
%
% MAIN PARAMETERS TO CHANGE
% fit_results_file:
%   Which Hopf fit result file to use.
%
% Cfg.TR:
%   Repetition time in seconds.
%
% Cfg.filt.lb / Cfg.filt.ub:
%   Bandpass filter cutoffs in Hz.
%
% Cfg.NSIM:
%   Number of baseline/perturbed simulations per subject-condition.
%
% Cfg.a / Cfg.sig:
%   Hopf bifurcation parameter and noise amplitude.
%
% Cfg.burn_in_seconds:
%   Burn-in duration before baseline simulation sampling.
%
% RUN ORDER
% Run after:
%   hopf_model_fit.m or a merged Hopf fit result file.
%
% Run before:
%   hopf_capacity_susceptibility_violin_plots.m

clear;
clc;
% LOAD FIT RESULTS

fit_results_file = 'hopf_sc_fit_results_schaefer1000_G0_to_4p5.mat';
if ~isfile(fit_results_file)
    error('Cannot find %s. Run hopf_model_fit.m first or change fit_results_file.', ...
        fit_results_file);
end
fit = load(fit_results_file);

OptimalG = fit.OptimalG;
Fdiff    = fit.Fdiff;
Tmax_all = fit.Tmax;

subjects_clean = fit.subjects_clean;
condition_order = fit.condition_order;

nSub  = size(OptimalG,1);
nCond = size(OptimalG,2);
% CFG

Cfg = struct();

Cfg.TR = 2.9;

Cfg.filt.lb = 0.008;
Cfg.filt.ub = 0.08;

Cfg.NSIM = 10;

Cfg.a = -0.02;
Cfg.sig = 0.01;

Cfg.burn_in_seconds = 2000;
Cfg.dt = 0.1 * Cfg.TR / 2;

rng(1);
% SC

if isfield(fit, 'SC')
    SC = double(fit.SC);
else
    SC = load_schaefer1000_sc();
end

SC(eye(size(SC))==1) = 0;

SC = SC ./ max(SC(:));
SC = SC * 0.2;

Cfg.nNodes = size(SC,1);

if Cfg.nNodes ~= 1000
    error('Expected a 1000 x 1000 SC matrix, found %d x %d.', size(SC,1), size(SC,2));
end

measure_weighting = 'SC';
% FILTER

fnq = 1/(2*Cfg.TR);

Wn = [Cfg.filt.lb/fnq Cfg.filt.ub/fnq];

[bfilt,afilt] = butter(2,Wn);
% OUTPUT

InfoCapacity  = NaN(nSub,nCond);
Susceptibility = NaN(nSub,nCond);
% MAIN LOOP

for cond = 1:nCond

    fprintf('\nCondition %d/%d\n',cond,nCond);

    for sub = 1:nSub

        fprintf(' Subject %d/%d\n',sub,nSub);

        G = OptimalG(sub,cond);

        f_diff = Fdiff{sub,cond};

        if numel(f_diff) ~= Cfg.nNodes
            error('Fdiff for subject %d condition %d has %d nodes, expected %d. Rerun hopf_model_fit.m at Schaefer-1000 first.', ...
                sub, cond, numel(f_diff), Cfg.nNodes);
        end

        Tmax = Tmax_all(sub,cond);

        ens_base = zeros(Cfg.NSIM,Cfg.nNodes);
        ens_pert = zeros(Cfg.NSIM,Cfg.nNodes);

        for isim = 1:Cfg.NSIM
            % BASELINE

            [ts,z_final] = simulate_hopf( ...
                SC,G,f_diff,Tmax,Cfg);

            ts = filter_bold(ts,bfilt,afilt);

            Ph = phase(ts);

            Ets = enstrophy_timeseries( ...
                Ph,SC);

            ens_base(isim,:) = ...
                mean(Ets,2,'omitnan')';
            % PERTURBED

            ts = simulate_hopf_pert( ...
                SC,G,f_diff,Tmax,Cfg,z_final);

            ts = filter_bold(ts,bfilt,afilt);

            Ph = phase(ts);

            Ets = enstrophy_timeseries( ...
                Ph,SC);

            ens_pert(isim,:) = ...
                mean(Ets,2,'omitnan')';

        end

       meanBaseline = mean(ens_base,1,'omitnan');

D = ens_pert - ...
    ones(Cfg.NSIM,1)*meanBaseline;

InfoCapacity(sub,cond) = ...
    mean(std(D,0,1,'omitnan'),'omitnan');

Susceptibility(sub,cond) = ...
    mean(mean(D,1,'omitnan'),'omitnan');

    end

end
% REPORT

fprintf('\n================ RESULTS ================\n');

for cond = 1:nCond

    fprintf('\n%s\n', string(condition_order{cond}));

    fprintf('Mean Info Capacity  = %.6f\n', ...
        mean(InfoCapacity(:,cond),'omitnan'));

    fprintf('Std  Info Capacity  = %.6f\n', ...
        std(InfoCapacity(:,cond),'omitnan'));

    fprintf('Mean Susceptibility = %.6f\n', ...
        mean(Susceptibility(:,cond),'omitnan'));

    fprintf('Std  Susceptibility = %.6f\n', ...
        std(Susceptibility(:,cond),'omitnan'));

end

fprintf('\n=========================================\n');
% SAVE

save('hopf_subjectwise_IC_SUS.mat', ...
    'InfoCapacity', ...
    'Susceptibility', ...
    'OptimalG', ...
    'subjects_clean', ...
    'condition_order', ...
    'measure_weighting', ...
    '-v7.3');

fprintf('DONE\n');
% HOPF SIM (UNCHANGED)

function [ts,z_final] = simulate_hopf(C,G,f_diff,Tmax,Cfg)

n = Cfg.nNodes;

dt = Cfg.dt;
dsig = sqrt(dt)*Cfg.sig;

omega = repmat(2*pi*f_diff(:),1,2);
omega(:,1) = -omega(:,1);

wC = G * C;
sumC = repmat(sum(wC,2),1,2);

a = Cfg.a * ones(n,2);
z = 0.1 * ones(n,2);

nBurn = round(Cfg.burn_in_seconds / dt);

for i = 1:nBurn
    z = hopf_step(z,a,omega,wC,sumC,dt,dsig);
end

step = round(Cfg.TR/dt);

ts = zeros(n,Tmax);

k = 0;

for i = 1:(Tmax*step)

    z = hopf_step( ...
        z,a,omega,wC,sumC,dt,dsig);

    if mod(i,step)==0

        k = k + 1;

        ts(:,k) = z(:,1);

    end

end

z_final = z;

end
% PERTURBED HOPF (UNCHANGED)

function ts = simulate_hopf_pert(C,G,f_diff,Tmax,Cfg,z0)

n = Cfg.nNodes;

dt = Cfg.dt;
dsig = sqrt(dt)*Cfg.sig;

omega = repmat(2*pi*f_diff(:),1,2);
omega(:,1) = -omega(:,1);

wC = G * C;
sumC = repmat(sum(wC,2),1,2);

a = -0.02 + 0.02*rand(n,1);
a = repmat(a,1,2);

z = z0;

step = round(Cfg.TR / dt);

ts = zeros(n,Tmax);
k = 0;

for i = 1:(Tmax*step)
    z = hopf_step(z,a,omega,wC,sumC,dt,dsig);
    if mod(i,step)==0
        k = k + 1;
        ts(:,k) = z(:,1);
    end
end

end
% HOPF STEP (UNCHANGED)

function z = hopf_step(z,a,omega,wC,sumC,dt,dsig)

suma = wC*z - sumC.*z;
zz = z(:,end:-1:1);

dz = a.*z + zz.*omega - z.*(z.^2 + zz.^2) + suma;

z = z + dt*dz + dsig*randn(size(z));

end
% FILTER (UNCHANGED)

function X = filter_bold(X,bfilt,afilt)

n = size(X,1);

for i = 1:n

    x = X(i,:) - mean(X(i,:));

    x = detrend(x);

    X(i,:) = filtfilt(bfilt,afilt,x);

end

end
% PHASE (UNCHANGED)

function Ph = phase(ts)

n = size(ts,1);
T = size(ts,2);

Ph = zeros(n,T);

for i = 1:n
    x = ts(i,:) - mean(ts(i,:));
    Ph(i,:) = angle(hilbert(x));
end

end
% SC-WEIGHTED ENSTROPHY

function E = enstrophy_timeseries(Ph,W)

n = size(Ph,1);
T = size(Ph,2);

z = complex(cos(Ph),sin(Ph));

E = zeros(n,T);

for i = 1:n

    w = W(i,:);
    w(i)=0;

    denom = sum(w);

    if denom==0
        E(i,:) = NaN;
        continue
    end

    local_order = (w * z) ./ denom;

    E(i,:) = abs(local_order);

end

end
