% hopf_model_fit.m
%
% PURPOSE
% Fit an SC-coupled Hopf model to empirical BOLD functional connectivity as
% a function of inter-parcel distance. For each subject, condition, and
% candidate global coupling G, the script simulates BOLD from the Hopf model,
% computes simulated FC-distance profiles, compares them with empirical
% FC-distance profiles, and stores the G with minimum fit error.
%
% INPUT FILES
% timeseries_data_all.mat
%   Created by create_timeseries_data_schaefer1000.py. Must contain Data,
%   subjects, and optionally condition_order.
%
% schaefer1000_7N_coords_decolab.csv
%   Schaefer-1000 parcel coordinates used for distance bins.
%
% SC_schaefer1000_decolab.mat or sc_schaefer1000_decolab.npy
%   Structural connectivity matrix loaded by load_schaefer1000_sc.m.
%
% REQUIRED HELPER FILES
% load_schaefer_coords.m
% load_schaefer1000_sc.m
% readNPY.m, only if loading SC from .npy
%
% OUTPUT FILES
% Cfg.save_file
%   Checkpoint file saved every Cfg.save_every_n_subjects.
%
% fit_results_file
%   Final Hopf fit result file containing Cfg, G_range, OptimalG,
%   OptimalError, FitErrorMean, FitErrorStd, empirical targets, Fdiff, Tmax,
%   SC, CoG, and optional final simulations.
%
% MAIN PARAMETERS TO CHANGE
% Cfg.filt.lb / Cfg.filt.ub:
%   Bandpass filter cutoffs in Hz.
%
% Cfg.NR:
%   Number of distance bins for FC-distance profiles.
%
% Cfg.fit_r_range_mm:
%   Distance-bin range used for fitting error.
%
% Cfg.TR:
%   Repetition time in seconds.
%
% Cfg.Glower / Cfg.Gstep / Cfg.Gupper:
%   Global coupling search grid.
%
% Cfg.NSIM_FIT:
%   Number of stochastic simulations per subject-condition-G.
%
% Cfg.NSIM_FINAL / Cfg.run_final_simulations:
%   Optional final simulations at each optimal G.
%
% Cfg.a:
%   Hopf bifurcation parameter.
%
% Cfg.sig:
%   Noise amplitude.
%
% RUN ORDER
% Run after:
%   create_timeseries_data_schaefer1000.py
%
% Run before:
%   hopf_model_fit_plots.m
%   hopf_capacity_susceptibility.m

Cfg.filt.lb = 0.008;
Cfg.filt.ub = 0.08;

Cfg.NR = 20;
Cfg.fit_r_range_mm = [8.13 33.82];

Cfg.nNodes = 1000; 
Cfg.TR = 2.9; 
Cfg.min_T = 10;

Cfg.Glower = 0;
Cfg.Gstep = 0.25;
Cfg.Gupper = 4.5;

Cfg.NSIM_FIT = 10;
Cfg.NSIM_FINAL = 100;
Cfg.run_final_simulations = false;

Cfg.a = -0.02;
Cfg.sig = 0.01;

Cfg.excTp = 0;

rng(1);

G_range = Cfg.Glower:Cfg.Gstep:Cfg.Gupper;
% AUTO-SAVE CONFIG

Cfg.save_file = 'hopf_checkpoint_schaefer1000_G0_to_4p5.mat';
Cfg.save_every_n_subjects = 1;   % you can set 1–5
fit_results_file = 'hopf_sc_fit_results_schaefer1000_G0_to_4p5.mat';
% LOAD DATA

input_file = 'timeseries_data_all.mat';
if ~isfile(input_file)
    error('Cannot find %s. Run create_timeseries_data_schaefer1000.py first.', ...
        input_file);
end

data = load(input_file);
Data = data.Data;
subjects = data.subjects;

if isfield(data, 'condition_order')
    condition_order = data.condition_order;
else
    condition_order = compose("Condition_%d", 1:size(Data,2));
end

nSub = size(Data,1);
nCond = size(Data,2);

keep_subject = true(nSub,1);
for sub = 1:nSub
    for cond = 1:nCond
        if isempty(Data{sub,cond}) || size(Data{sub,cond},1) < Cfg.min_T
            keep_subject(sub) = false;
            break;
        end
    end
end

Data_clean = Data(keep_subject,:);
subjects_clean = subjects(keep_subject);

nSub_clean = size(Data_clean,1);
nCond_clean = size(Data_clean,2);

fprintf('Kept subjects = %d / %d\n', nSub_clean, nSub);
fprintf('Conditions = %d\n', nCond_clean);
% LOAD COORDINATES AND SC

[CoG, ~] = load_schaefer_coords('schaefer1000_7N_coords_decolab.csv');

SC = load_schaefer1000_sc();

if size(SC,1) ~= Cfg.nNodes || size(SC,2) ~= Cfg.nNodes
    error('SC matrix must be %d x %d, found %d x %d.', ...
        Cfg.nNodes, Cfg.nNodes, size(SC,1), size(SC,2));
end

if size(CoG,1) ~= Cfg.nNodes
    error('Coordinate file must contain %d parcels, found %d.', ...
        Cfg.nNodes, size(CoG,1));
end

SC(eye(Cfg.nNodes) == 1) = 0;
SC = SC ./ max(SC(:));
SC = SC * 0.2;

Cfg.nNodes = size(SC,1);

fprintf('Using Schaefer-1000 SC: %d x %d\n', size(SC,1), size(SC,2));
fprintf('SC max after normalization = %.4f\n', max(SC(:)));

if size(Data_clean{1,1}, 2) ~= Cfg.nNodes
    error(['Expected time-series data with %d parcels, found %d. ' ...
           'Replace timeseries_data_all.mat with Schaefer-1000 BOLD data before running hopf_model_fit.m.'], ...
           Cfg.nNodes, size(Data_clean{1,1}, 2));
end
% DISTANCE BINS

rr = squareform(pdist(CoG));
distance_range = max(rr(:));
delta = distance_range / Cfg.NR;
distance_bin_edges = linspace(0, distance_range, Cfg.NR + 1);
distance_bin_centers = distance_bin_edges(1:end-1) + diff(distance_bin_edges)/2;
fitBinMask = distance_bin_centers >= Cfg.fit_r_range_mm(1) & ...
    distance_bin_centers <= Cfg.fit_r_range_mm(2);

if ~any(fitBinMask)
    error('No distance bins found inside requested fitting range %.2f-%.2f mm.', ...
        Cfg.fit_r_range_mm(1), Cfg.fit_r_range_mm(2));
end

Cfg.NRini = find(fitBinMask, 1, 'first');
Cfg.NRfin = find(fitBinMask, 1, 'last');
Cfg.fit_bin_range_mm = [ ...
    distance_bin_edges(Cfg.NRini), ...
    distance_bin_edges(Cfg.NRfin + 1)];

fprintf('Requested fitting range r = %.2f-%.2f mm\n', ...
    Cfg.fit_r_range_mm(1), Cfg.fit_r_range_mm(2));
fprintf('Using distance bins %d-%d, covering %.2f-%.2f mm\n', ...
    Cfg.NRini, Cfg.NRfin, ...
    Cfg.fit_bin_range_mm(1), Cfg.fit_bin_range_mm(2));
% FILTER

fnq = 1 / (2 * Cfg.TR);
Wn = [Cfg.filt.lb / fnq, Cfg.filt.ub / fnq];
[bfilt, afilt] = butter(2, Wn);
% EMPIRICAL TARGETS

EmpCorrFcn = cell(nSub_clean, nCond_clean);
EmpFC = cell(nSub_clean, nCond_clean);
Fdiff = cell(nSub_clean, nCond_clean);
Tmax = NaN(nSub_clean, nCond_clean);

fprintf('\nPreparing empirical FC-distance targets...\n');

for sub = 1:nSub_clean
    for cond = 1:nCond_clean

        BOLD = Data_clean{sub,cond}';
        BOLD_filt = filter_bold(BOLD, bfilt, afilt, Cfg.excTp);

        FC = corrcoef(BOLD_filt');
        EmpFC{sub,cond} = FC;

        EmpCorrFcn{sub,cond} = fc_distance_profile( ...
            FC, rr, Cfg.NR, delta);

        Fdiff{sub,cond} = estimate_node_frequencies_ttb( ...
    BOLD_filt, Cfg.TR);

        Tmax(sub,cond) = size(BOLD_filt,2);
    end
end
% FIT G

FitErrorMean = NaN(nSub_clean, nCond_clean, numel(G_range));
FitErrorStd  = NaN(nSub_clean, nCond_clean, numel(G_range));

OptimalG = NaN(nSub_clean, nCond_clean);
OptimalError = NaN(nSub_clean, nCond_clean);

fprintf('\n========================================\n');
fprintf('HOPF SC FIT STARTED\n');
fprintf('Bifurcation parameter a = %.3f\n', Cfg.a);
fprintf('Simulations per subject-condition-G = %d\n', Cfg.NSIM_FIT);
fprintf('G range = %.3f:%.3f:%.3f\n', ...
    Cfg.Glower, Cfg.Gstep, Cfg.Gupper);
fprintf('========================================\n');

for cond = 1:nCond_clean

    fprintf('\n########################################\n');
    fprintf('CONDITION %d / %d: %s\n', ...
        cond, nCond_clean, string(condition_order{cond}));
    fprintf('########################################\n');

    for sub = 1:nSub_clean

        fprintf('\nSubject %d / %d: %s\n', ...
            sub, nSub_clean, string(subjects_clean{sub}));

        empcorrfcn = EmpCorrFcn{sub,cond};
        f_diff = Fdiff{sub,cond};
        thisTmax = Tmax(sub,cond);

        for iG = 1:numel(G_range)

            G = G_range(iG);
            sim_errors = NaN(Cfg.NSIM_FIT,1);

            for isim = 1:Cfg.NSIM_FIT

                sim_ts = simulate_hopf_sc( ...
                    SC, G, f_diff, thisTmax, Cfg);

                if any(~isfinite(sim_ts(:)))
                    fprintf(['    Skipping simulation %d: non-finite ' ...
                        'simulated BOLD values at G = %.3f | sub = %d | cond = %d\n'], ...
                        isim, G, sub, cond);
                    continue
                end

                sim_filt = filter_bold(sim_ts, bfilt, afilt, 0);

                simFC = corrcoef(sim_filt');

                simcorrfcn = fc_distance_profile( ...
                    simFC, rr, Cfg.NR, delta);

                sim_errors(isim) = profile_error( ...
                    simcorrfcn, empcorrfcn, ...
                    Cfg.NRini, Cfg.NRfin);
            end

            FitErrorMean(sub,cond,iG) = mean(sim_errors,'omitnan');
            FitErrorStd(sub,cond,iG)  = std(sim_errors,'omitnan');

            fprintf('  G = %.3f | error = %.6f +/- %.6f\n', ...
                G, ...
                FitErrorMean(sub,cond,iG), ...
                FitErrorStd(sub,cond,iG));
        end

        subjectErrors = squeeze(FitErrorMean(sub,cond,:));

        if any(isfinite(subjectErrors))
            [OptimalError(sub,cond), bestIdx] = ...
                min(subjectErrors, [], 'omitnan');

            OptimalG(sub,cond) = G_range(bestIdx);
        else
            OptimalError(sub,cond) = NaN;
            OptimalG(sub,cond) = NaN;

            warning(['No finite Hopf simulations for subject %s condition %s. ' ...
                'Optimal G left as NaN.'], ...
                string(subjects_clean{sub}), string(condition_order{cond}));
        end
        % AUTO-SAVE (NEW ADDITION)

        if mod(sub, Cfg.save_every_n_subjects) == 0
            fprintf('\n[SAVING CHECKPOINT]\n');

            save(Cfg.save_file, ...
                'FitErrorMean','FitErrorStd', ...
                'OptimalG','OptimalError', ...
                'G_range','Cfg', ...
                'subjects_clean','condition_order', ...
                '-v7.3');
        end


        fprintf(['  >>> OPTIMAL G subject %s condition %s = ', ...
                 '%.3f | error = %.6f\n'], ...
                 string(subjects_clean{sub}), ...
                 string(condition_order{cond}), ...
                 OptimalG(sub,cond), ...
                 OptimalError(sub,cond));
    end
end

MeanG_by_condition = mean(OptimalG,1,'omitnan');
StdG_by_condition  = std(OptimalG,0,1,'omitnan');
MeanG_all = mean(OptimalG(:),'omitnan');

fprintf('\n========================================\n');
fprintf('OPTIMAL G BY SUBJECT AND CONDITION\n');
fprintf('========================================\n');

disp(array2table( ...
    OptimalG, ...
    'VariableNames', matlab.lang.makeValidName(string(condition_order))));

fprintf('\n========================================\n');
fprintf('AVERAGE OPTIMAL G BY CONDITION\n');
fprintf('========================================\n');

for cond = 1:nCond_clean
    fprintf('%s | mean G = %.4f | std G = %.4f\n', ...
        string(condition_order{cond}), ...
        MeanG_by_condition(cond), ...
        StdG_by_condition(cond));
end

fprintf('Overall mean optimal G = %.4f\n', MeanG_all);
% OPTIONAL FINAL SIMULATIONS

FinalSimFC = cell(nSub_clean, nCond_clean);
FinalSimCorrFcn = cell(nSub_clean, nCond_clean);

if Cfg.run_final_simulations

    fprintf('\nRunning final simulations at each subject-condition optimal G...\n');

    for cond = 1:nCond_clean
        for sub = 1:nSub_clean

            G = OptimalG(sub,cond);
            f_diff = Fdiff{sub,cond};
            thisTmax = Tmax(sub,cond);

            if ~isfinite(G)
                fprintf(['    Skipping final simulations: no finite optimal G ' ...
                    'for sub = %d | cond = %d\n'], sub, cond);
                continue
            end

            fc_stack = NaN(Cfg.nNodes, Cfg.nNodes, Cfg.NSIM_FINAL);
            corr_stack = NaN(Cfg.nNodes, Cfg.NR, Cfg.NSIM_FINAL);

            for isim = 1:Cfg.NSIM_FINAL

                sim_ts = simulate_hopf_sc( ...
                    SC, G, f_diff, thisTmax, Cfg);

                if any(~isfinite(sim_ts(:)))
                    fprintf(['    Skipping final simulation %d: non-finite ' ...
                        'simulated BOLD values at G = %.3f | sub = %d | cond = %d\n'], ...
                        isim, G, sub, cond);
                    continue
                end

                sim_filt = filter_bold(sim_ts, bfilt, afilt, 0);

                simFC = corrcoef(sim_filt');

                fc_stack(:,:,isim) = simFC;

                corr_stack(:,:,isim) = fc_distance_profile( ...
                    simFC, rr, Cfg.NR, delta);
            end

            FinalSimFC{sub,cond} = mean(fc_stack,3,'omitnan');
            FinalSimCorrFcn{sub,cond} = mean(corr_stack,3,'omitnan');
        end
    end
end

save(fit_results_file, ...
    'Cfg', ...
    'G_range', ...
    'OptimalG', ...
    'OptimalError', ...
    'MeanG_by_condition', ...
    'StdG_by_condition', ...
    'MeanG_all', ...
    'FitErrorMean', ...
    'FitErrorStd', ...
    'subjects_clean', ...
    'condition_order', ...
    'EmpCorrFcn', ...
    'EmpFC', ...
    'Fdiff', ...
    'Tmax', ...
    'SC', ...
    'CoG', ...
    'FinalSimFC', ...
    'FinalSimCorrFcn', ...
    '-v7.3');

fprintf('\nSaved: %s\n', fit_results_file);
% LOCAL FUNCTIONS

function BOLD_filt = filter_bold(BOLD, bfilt, afilt, excTp)

    nAreas = size(BOLD,1);
    nTime = size(BOLD,2);

    BOLD_filt = zeros(nAreas, nTime - 2*excTp);

    for seed = 1:nAreas

        sig = BOLD(seed,:) - mean(BOLD(seed,:));
        sig = detrend(sig);
        sig = filtfilt(bfilt, afilt, sig);

        if excTp > 0
            sig = sig(excTp+1:end-excTp);
        end

        BOLD_filt(seed,:) = sig;
    end
end

function f_diff = estimate_node_frequencies_ttb(BOLD_filt,TR)

[nAreas,TT] = size(BOLD_filt);

freq = (0:floor(TT/2)-1)/(TT*TR);

PowSpect = zeros(length(freq),nAreas);

for seed = 1:nAreas

    pw = abs(fft(BOLD_filt(seed,:)));

    PowSpect(:,seed) = ...
        pw(1:floor(TT/2)).^2/(TT/TR);

end

for seed = 1:nAreas

    PowSpect(:,seed) = ...
        gaussfilt(freq, PowSpect(:,seed)',0.01);

end

[~,idx] = max(PowSpect);

f_diff = freq(idx)';

f_diff(f_diff==0) = ...
    mean(f_diff(f_diff~=0));

end

function corrfcn = fc_distance_profile(FC, rr, NR, delta)

    nAreas = size(FC,1);
    corrfcn = NaN(nAreas, NR);

    for i = 1:nAreas

        numind = zeros(1,NR);
        corrfcn_1 = zeros(1,NR);

        for j = 1:nAreas

            r = rr(i,j);
            index = floor(r / delta) + 1;

            if index == NR + 1
                index = NR;
            end

            mcc = FC(i,j);

            if ~isnan(mcc)
                corrfcn_1(index) = corrfcn_1(index) + mcc;
                numind(index) = numind(index) + 1;
            end
        end

        corrfcn(i,:) = corrfcn_1 ./ numind;
    end
end

function err = profile_error( ...
    simcorrfcn,empcorrfcn,NRini,NRfin)

    nNodes = size(simcorrfcn,1);

    err1 = NaN(nNodes,1);

    for i = 1:nNodes

        err11 = NaN(1,size(simcorrfcn,2));

        for k = NRini:NRfin

            err11(k) = ...
                (simcorrfcn(i,k)-empcorrfcn(i,k))^2;

        end

        err1(i) = ...
            mean(err11(NRini:NRfin),'omitnan');

    end

    err = sqrt(mean(err1,'omitnan'));

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

    %% Burn-in exactly as TTB

    for t = 0:dt:2000

        z = hopf_step( ...
            z,a,omega,wC,sumC,dt,dsig);

    end

    %% Main simulation

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

function z = hopf_step(z, a, omega, wC, sumC, dt, dsig)

    suma = wC*z - sumC.*z;
    zz = z(:,end:-1:1);

    z = z + dt * ( ...
        a.*z + ...
        zz.*omega - ...
        z.*(z.*z + zz.*zz) + ...
        suma ) + ...
        dsig * randn(size(z));
end

function y = gaussfilt(x,yin,sigma)

    x = x(:);
    yin = yin(:);

    n = length(x);
    y = zeros(n,1);

    for i = 1:n

        w = exp(-(x - x(i)).^2/(2*sigma^2));
        w = w/sum(w);

        y(i) = sum(w .* yin);

    end
end
