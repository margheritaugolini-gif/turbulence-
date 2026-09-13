% BOLD2hilbertupdated.m
%
% PURPOSE
% 
% Main preprocessing script for the Schaefer-1000 workflow. It loads the
% subject/condition BOLD time series, removes subjects with missing or very
% short recordings, bandpass-filters each parcel signal, computes Hilbert
% phases, computes static functional connectivity, and saves outputs used by
% the later Kuramoto/turbulence analyses.
%
% INPUT FILE
% 
% timeseries_data_all.mat
%   Created by create_timeseries_data_schaefer1000.py.
%   Must contain:
%     Data     : cell array, subjects x conditions, each cell time x parcels
%     subjects : subject labels
%
% OUTPUT FILES
% 
% preprocessed_data.mat
%   Saves FC, Phase, filtered BOLD_out, subjects_clean, Tmap, TBL, and
%   preprocessing parameters.
%
% phase_for_turbulence.mat
%   Smaller file used by whole_brain_turbulence.m and network_node_turbulence.m.
%   Saves Phase_all plus subject/time/audit metadata.
%
% MAIN PARAMETERS TO CHANGE
% 
% min_T    : minimum accepted number of time points per subject-condition.
% TR       : repetition time in seconds.
% low_cut  : lower bandpass cutoff in Hz.
% high_cut : upper bandpass cutoff in Hz.
% excTp    : number of time points removed from each edge after filtering.
%
% RUN ORDER
% 
% Run after:
%   create_timeseries_data_schaefer1000.py
%
% Run before:
%   metastability_kuramoto.m
%   whole_brain_turbulence.m
%   network_node_turbulence.m
%
% NOTES
%
% The script expects 1000 parcels. It stops with an error if the input data
% are not Schaefer-1000 time series.


% LOAD DATA
clear; clc;

input_file = 'timeseries_data_all.mat';
if ~isfile(input_file)
    error('Cannot find %s. Run create_timeseries_data_schaefer1000.py first.', ...
        input_file);
end

data = load(input_file);

Data = data.Data;
subjects = data.subjects;

nSub  = size(Data,1);
nCond = size(Data,2);

min_T = 10;

% SUBJECT EXCLUSION 
keep_subject = true(nSub,1);
exclusion_reason = strings(nSub,1);

for i = 1:nSub
    
    for j = 1:nCond
        
        % missing data check
        if isempty(Data{i,j})
            keep_subject(i) = false;
            exclusion_reason(i) = "missing cell at condition " + j;
            break;
        end
        
        % minimum time length check
        if size(Data{i,j},1) < min_T
            keep_subject(i) = false;
            exclusion_reason(i) = "T < min_T at condition " + j;
            break;
        end
        
    end
    
    if keep_subject(i)
        exclusion_reason(i) = "kept";
    end
    
end


% APPLY CLEANING
Data_clean = Data(keep_subject,:);
subjects_clean = subjects(keep_subject);

removed_subjects = subjects(~keep_subject);

fprintf('Final N = %d\n', sum(keep_subject));

disp('Removed subjects:')
disp(removed_subjects)

SubjectID = (1:nSub)';

TBL = table(SubjectID, keep_subject, exclusion_reason);

disp('=== SUBJECT EXCLUSION REPORT ===')
disp(TBL)


% T-MAP CONSTRUCTION

nSub_clean  = size(Data_clean,1);
nCond_clean = size(Data_clean,2);

Tmap = NaN(nSub_clean, nCond_clean);

for i = 1:nSub_clean
    for j = 1:nCond_clean
        Tmap(i,j) = size(Data_clean{i,j},1);
    end
end

disp('Tmap (subjects × conditions):')
disp(Tmap);

fprintf('\nT stats:\n');
fprintf('min T = %d\n', min(Tmap(:)));
fprintf('max T = %d\n', max(Tmap(:)));
fprintf('mean T = %.2f\n', mean(Tmap(:)));


% FILTER 

TR = 2.9;
fs = 1/TR;

low_cut  = 0.008;
high_cut = 0.08;

[bfilt, afilt] = butter(2, [low_cut high_cut]/(fs/2));

excTp = 3;
nAreas = size(Data_clean{1,1}, 2);

if nAreas ~= 1000
    error(['Expected 1000 regions in timeseries_data_all.mat, found %d. ' ...
           'Replace the time-series input with Schaefer-1000 data before preprocessing.'], ...
           nAreas);
end



% INITIALIZE OUTPUTS

FC = cell(nSub_clean, nCond_clean);
Phase = cell(nSub_clean, nCond_clean);
BOLD_out = cell(nSub_clean, nCond_clean);



% PREPROCESSING 

for i = 1:nSub_clean
    for j = 1:nCond_clean
        
        % regions x T format
        BOLD = Data_clean{i,j}';

        if size(BOLD,1) ~= nAreas
            error('Expected %d regions, found %d for subject %d condition %d.', ...
                nAreas, size(BOLD,1), i, j);
        end
        
        nTime = size(BOLD,2);
        
        BOLD_processed = zeros(nAreas, nTime - 2*excTp);
        Phase_BOLD = zeros(nAreas, nTime - 2*excTp);
        
        for seed = 1:nAreas
            
            sig = BOLD(seed,:) - mean(BOLD(seed,:));
            sig = filtfilt(bfilt, afilt, sig);
            sig = sig(excTp+1:end-excTp);
            
            BOLD_processed(seed,:) = sig;
            Phase_BOLD(seed,:) = angle(hilbert(sig));
            
        end
        
        % static FC
        FC{i,j} = corrcoef(BOLD_processed');
        Phase{i,j} = Phase_BOLD;
        BOLD_out{i,j} = BOLD_processed;
        
    end
end


% SAVE OUTPUT

Phase_all = Phase;

save('preprocessed_data.mat', ...
    'FC', 'Phase', 'BOLD_out', ...
    'subjects_clean', 'Tmap', 'TBL', ...
    'low_cut', 'high_cut', 'TR', 'excTp');

save('phase_for_turbulence.mat', ...
    'Phase_all', 'subjects_clean', 'Tmap', 'TBL', ...
    'low_cut', 'high_cut', 'TR', 'excTp');
