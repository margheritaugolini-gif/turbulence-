% network_node_turbulence_helper.m
%
% PURPOSE
% Helper function for node-level turbulence. Given Hilbert phases and parcel
% coordinates, it computes one turbulence value per parcel, subject, and
% condition using local Kuramoto dynamics at a fixed spatial scale.
%
% INPUTS
% Phase_all
%   Cell array, subjects x conditions. Each cell is regions x time and
%   contains Hilbert phases.
%
% CoG
%   Coordinate matrix, regions x 3.
%
% OUTPUT
% output.NodeTurbulence
%   Cell array, subjects x conditions. Each cell is regions x 1.
%
% output.R
%   Local Kuramoto order-parameter values.
%
% output.lambda
%   Fixed spatial scale used for node turbulence.
%
% MAIN PARAMETERS TO CHANGE
% lambda:
%   Fixed spatial scale. Current setting is 0.12.
%
% USED BY
% network_node_turbulence.m



function output = network_node_turbulence_helper(Phase_all, CoG)
% BASIC INFORMATION

NPARCELS = size(CoG,1);

NSUB = size(Phase_all,1);

NCOND = size(Phase_all,2);

lambda = 0.12;

rr = zeros(NPARCELS,NPARCELS);


for i = 1:NPARCELS
    
    for j = 1:NPARCELS
        
        rr(i,j) = norm(CoG(i,:) - CoG(j,:));
        
    end
end

% CREATE SPATIAL WEIGHT MATRIX

C = exp(-lambda * rr);


NodeTurbulence = cell(NSUB,NCOND);

R_all = cell(NSUB,NCOND);


% MAIN ANALYSIS LOOP

for sub = 1:NSUB
    
    fprintf('\n');
    fprintf('====================================\n');
    fprintf('SUBJECT %d / %d\n',sub,NSUB);
    fprintf('====================================\n');

    % LOOP ACROSS CONDITIONS
    
    for cond = 1:NCOND
        
        fprintf('Condition %d / %d\n',cond,NCOND);
        % SKIP EMPTY CELLS
        
        if isempty(Phase_all{sub,cond})
            
            fprintf('Condition empty -> skipped\n');
            
            continue
        end

        % GET PHASE DATA
        
        Phases = Phase_all{sub,cond};
    
        
        Tmax = size(Phases,2);
        
        % Number of timepoints
        % CONVERT PHASES TO COMPLEX FORM
        
        complex_phase = exp(1i * Phases);
       

        % LOCAL KURAMOTO ORDER PARAMETER
        
        R = zeros(NPARCELS,Tmax);
        
        % R(region,time)
        %
        % measures local synchronization
        
        
        for i = 1:NPARCELS
            % GET SPATIAL WEIGHTS
            
            weights = C(i,:)';
            
            % [200 x 1]
            % APPLY WEIGHTS TO PHASES
            
            weighted_phase = ...
                repmat(weights,1,Tmax) .* complex_phase;
            % COMPUTE WEIGHTED MEAN PHASE
            
            sumphases = ...
                sum(weighted_phase,1) / sum(weights);
            % LOCAL SYNCHRONIZATION
            
            R(i,:) = abs(sumphases);
            
           
            
        end
        % NODE-LEVEL TURBULENCE
        
        NodeTurb = zeros(NPARCELS,1);
        
        
        for i = 1:NPARCELS
            
            NodeTurb(i) = std(R(i,:));
           
            
        end
      
        
        NodeTurbulence{sub,cond} = NodeTurb;
        
        R_all{sub,cond} = R;
        
        
        fprintf('Finished condition %d\n',cond);
        
    end
end

% FINAL OUTPUT

output.NodeTurbulence = NodeTurbulence;

output.R = R_all;

output.lambda = lambda;



fprintf('\n');
fprintf('====================================\n');
fprintf('NODE TURBULENCE COMPLETED\n');
fprintf('====================================\n');

end
