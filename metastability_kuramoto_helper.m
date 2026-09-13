% metastability_kuramoto_helper.m
%
% PURPOSE
% -------
% Helper function for global Kuramoto summaries. Computes the Kuramoto
% order-parameter time series, mean synchronization, and metastability from
% a regions-by-time Hilbert phase matrix.

function [OP, Synchro, Metasta] = metastability_kuramoto_helper(Phase_BOLD, nAreas)

T = size(Phase_BOLD,2); %How many time points do we have, counting basically 

OP = zeros(1,T); %Prepare empty storage for synchronization over time (kuramoto)

for t = 1:T
    OP(t) = abs(mean(exp(1i * Phase_BOLD(:,t)))); %for each TIME POINT, compute synchronization across all brain regions
end

Synchro = mean(OP); %Overall synchronization of the brain for time points
Metasta = std(OP); %overall metastability 

end
