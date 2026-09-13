function SC = load_schaefer1000_sc()
%LOAD_SCHAEFER1000_SC Load the Schaefer-1000 structural connectome.

matFile = 'SC_schaefer1000_decolab.mat';
npyFile = 'sc_schaefer1000_decolab.npy';

if exist(matFile, 'file')
    scData = load(matFile, 'SC');
    SC = scData.SC;
elseif exist(npyFile, 'file')
    SC = readNPY(npyFile);
else
    error('Cannot find %s or %s.', matFile, npyFile);
end

SC = double(SC);
end
