% whole_brain_turbulence.m
%
% PURPOSE
% -------
% Run the whole-brain turbulence computation across all subjects and
% conditions. This script loads Hilbert phases and Schaefer-1000 parcel
% coordinates, calls whole_brain_turbulence_helper.m, and saves the full
% turbulence output.
%
% INPUT FILES
% -----------
% phase_for_turbulence.mat
%   Created by BOLD2hilbertupdated.m. Must contain Phase_all.
%
% schaefer1000_7N_coords_decolab.csv
%   Parcel coordinate file. The number of rows must match the number of
%   regions in Phase_all.
%
% REQUIRED HELPER FILES
% ---------------------
% whole_brain_turbulence_helper.m
% load_schaefer_coords.m
%
% OUTPUT FILE
% -----------
% turbulence_results.mat
%   Saves output.Turbulence, output.R, and output.LAMBDA.
%
% MAIN PARAMETERS TO CHANGE
% -------------------------
% Coordinate filename:
%   passed to load_schaefer_coords(...)
%
% Lambda range:
%   set inside whole_brain_turbulence_helper.m.
%
% RUN ORDER
% ---------
% Run after:
%   BOLD2hilbertupdated.m
%
% Run before:
%   whole_brain_turbulence_violin_plots.m
%   information_cascade.m

clear;
clc;

% SCRIPT TO RUN TURBULENCE ANALYSIS


% LOAD PHASE DATA

input_file = 'phase_for_turbulence.mat';
if ~isfile(input_file)
    error('Cannot find %s. Run BOLD2hilbertupdated.m first.', input_file);
end

load(input_file);

% Must contain:
%
% Phase_all
%
% Structure:
% {subject , condition}
%
% Each cell:
% [regions x time]



% LOAD BRAIN REGION COORDINATES

[CoG, ~] = load_schaefer_coords('schaefer1000_7N_coords_decolab.csv');



% EXTRACT XYZ COORDINATES

% Creates:
%
% [regions x 3]
%
% Columns:
% X Y Z coordinates



% CHECK DATA MATCHES

if size(CoG,1) ~= size(Phase_all{1,1},1)
    
    error(['Mismatch between number of brain regions ' ...
           'in coordinates and phase data']);
end



% RUN TURBULENCE COMPUTATION

output = whole_brain_turbulence_helper(Phase_all, CoG);



% SAVE RESULTS

save('turbulence_results.mat','output','-v7.3');

disp('Turbulence analysis completed.');
