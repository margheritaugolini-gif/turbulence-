% network_node_turbulence.m
%
% PURPOSE
% Run node-level turbulence for every subject and condition. This script
% loads Hilbert phases and Schaefer-1000 parcel coordinates, calls
% network_node_turbulence_helper.m, and saves one turbulence value per
% parcel, subject, and condition.
%
% INPUT FILES
% phase_for_turbulence.mat
%   Created by BOLD2hilbertupdated.m. Must contain Phase_all.
%
% schaefer1000_7N_coords_decolab.csv
%   Parcel coordinate file. The number of rows must match the number of
%   regions in Phase_all.
%
% REQUIRED HELPER FILES
% network_node_turbulence_helper.m
% load_schaefer_coords.m
%
% OUTPUT FILE
% node_turbulence_results.mat
%   Saves output.NodeTurbulence, output.R, and output.lambda.
%
% MAIN PARAMETERS TO CHANGE
% Coordinate filename:
%   passed to load_schaefer_coords(...)
%
% Fixed lambda:
%   set inside network_node_turbulence_helper.m.
%
% RUN ORDER
% Run after:
%   BOLD2hilbertupdated.m
%
% Run before:
%   network_turbulence_spider_plots.m, if using node_turbulence_results.mat.

clear;
clc;

% LOAD PHASE DATA

input_file = 'phase_for_turbulence.mat';
if ~isfile(input_file)
    error('Cannot find %s. Run BOLD2hilbertupdated.m first.', input_file);
end

load(input_file);

% LOAD SCHAEFER COORDINATES

[CoG, ~] = load_schaefer_coords('schaefer1000_7N_coords_decolab.csv');

% CHECK DATA MATCHES

if size(CoG,1) ~= size(Phase_all{1,1},1)
    
    error(['Mismatch between number of regions ' ...
           'in coordinates and phase data']);
end
% DISPLAY BASIC INFORMATION

NSUB = size(Phase_all,1);

NCOND = size(Phase_all,2);

NPARCELS = size(CoG,1);


fprintf('\n');
fprintf('====================================\n');
fprintf('DATA INFORMATION\n');
fprintf('====================================\n');

fprintf('Subjects   : %d\n',NSUB);

fprintf('Conditions : %d\n',NCOND);

fprintf('Parcels    : %d\n',NPARCELS);

% RUN NODE-LEVEL TURBULENCE

fprintf('\n');
fprintf('====================================\n');
fprintf('RUNNING NODE TURBULENCE\n');
fprintf('====================================\n');


output = network_node_turbulence_helper(Phase_all, CoG);

% SAVE RESULTS

save('node_turbulence_results.mat','output');

% FINISHED

fprintf('\n');
fprintf('====================================\n');
fprintf('RESULTS SAVED\n');
fprintf('====================================\n');

fprintf('File saved:\n');

fprintf('node_turbulence_results.mat\n');


disp('DONE');
