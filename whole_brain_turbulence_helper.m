% whole_brain_turbulence_helper.m
%
% PURPOSE
% -------
% Helper function for whole-brain turbulence. Given Hilbert phases and parcel
% coordinates, it computes local Kuramoto synchronization and turbulence
% across spatial lambda scales for every subject and condition.
%
% INPUTS
% ------
% Phase_all
%   Cell array, subjects x conditions. Each cell is regions x time and
%   contains Hilbert phases.
%
% CoG
%   Coordinate matrix, regions x 3.
%
% OUTPUT
% ------
% output.Turbulence
%   Cell array, subjects x conditions. Each cell is nLambda x 1.
%
% output.R
%   Local synchronization values, stored as lambda x regions x time.
%
% output.LAMBDA
%   Spatial scale values used for the analysis.
%
% MAIN PARAMETERS TO CHANGE
% -------------------------
% LAMBDA:
%   Spatial scale range. Current setting is 0.01:0.01:0.30.
%
% USED BY
% -------
% whole_brain_turbulence.m

function output = whole_brain_turbulence_helper(Phase_all, CoG)



% BASIC INFORMATION

NPARCELLS = size(CoG,1);
% Number of brain regions
% Example: 200 regions

NSUB = size(Phase_all,1);
% Number of subjects

NCOND = size(Phase_all,2);
% Number of experimental conditions



% DEFINE SPATIAL SCALES (LAMBDA)

LAMBDA = 0.01:0.01:0.30;

% Creates 30 values between:
% 0.01 → very broad/global interactions
% 0.30 → very local interactions

NLAMBDA = length(LAMBDA);



% COMPUTE DISTANCE BETWEEN ALL REGIONS

rr = zeros(NPARCELLS,NPARCELLS);

% rr(i,j) will contain:
% distance between region i and region j

for i = 1:NPARCELLS
    
    for j = 1:NPARCELLS
        
        % Euclidean distance between coordinates
        
        rr(i,j) = norm(CoG(i,:) - CoG(j,:));
        
    end
end



% CREATE SPATIAL WEIGHT MATRICES

C1 = zeros(NLAMBDA,NPARCELLS,NPARCELLS); %Each slice corresponds to a different decay parameter λ

% For each lambda:
% create a matrix of spatial weights

for ilam = 1:NLAMBDA
    
    lambda = LAMBDA(ilam);
    
    % Exponential decay:
    %
    % close regions → strong weight
    % far regions   → weak weight
    
    C1(ilam,:,:) = exp(-lambda * rr);
    
end



% STORAGE VARIABLES

Turbulence = cell(NSUB,NCOND);

% Will store:
% turbulence across lambda values

R_all = cell(NSUB,NCOND);

% Will store:
% local synchronization values



% MAIN ANALYSIS LOOP

for sub = 1:NSUB
    
    fprintf('Processing Subject %d / %d\n',sub,NSUB);
    
    
    for cond = 1:NCOND
        
        % Skip empty cells
        
        if isempty(Phase_all{sub,cond})
            continue
        end
        
        
        %% ------------------------------------------------
        % GET PHASE DATA
        % -------------------------------------------------
        
        Phases = Phase_all{sub,cond};
        
        % Size:
        % [regions x time]
        
        Tmax = size(Phases,2);
        
        
        %% ------------------------------------------------
        % CONVERT PHASES TO COMPLEX FORM
        % -------------------------------------------------
        
        complex_phase = exp(1i * Phases);
        
        % Converts:
        %
        % theta  → e^(i*theta)
        %
        % This places each phase on the unit circle
        %
        % Necessary for Kuramoto synchronization
        
        
        
        %% ------------------------------------------------
        % LOCAL KURAMOTO ORDER PARAMETER
        % -------------------------------------------------
        
        R_lambda = zeros(NLAMBDA,NPARCELLS,Tmax);
        
        % Dimensions:
        %
        % lambda x region x time
        
        
        for ilam = 1:NLAMBDA
            
            % Get weight matrix for this lambda
            
            C1lam = squeeze(C1(ilam,:,:));
            
            
            for i = 1:NPARCELLS
                
                % Weights from region i to all regions
                
                weights = C1lam(i,:)';
                
                
                % Multiply each region phase
                % by its spatial weight
                
                weighted_phase = ...
                    repmat(weights,1,Tmax) .* complex_phase;
                
                
                % Average weighted phases
                
                sumphases = ...
                    sum(weighted_phase,1) / sum(weights);
                
                
                % Magnitude = synchronization strength
                
                R_lambda(ilam,i,:) = abs(sumphases);
                
                % Values near:
                %
                % 1 → highly synchronized
                % 0 → desynchronized
                
            end
        end
        
        
        
        %% ------------------------------------------------
        % COMPUTE TURBULENCE
        % -------------------------------------------------
        
        T = zeros(NLAMBDA,1);
        
        
        for ilam = 1:NLAMBDA
            
            temp = squeeze(R_lambda(ilam,:,:));
            
            % Compute variability of synchronization
            
            T(ilam) = std(temp(:));
            
            % Higher std =
            % more fluctuations =
            % more turbulence
            
        end
        
        
        
        %% ------------------------------------------------
        % STORE RESULTS
        % -------------------------------------------------
        
        Turbulence{sub,cond} = T;
        
        R_all{sub,cond} = R_lambda;
        
    end
end



% FINAL OUTPUT

output.Turbulence = Turbulence;

output.R = R_all;

output.LAMBDA = LAMBDA;

end
