function [CoG, Labels] = load_schaefer_coords(filename)
%LOAD_SCHAEFER_COORDS Load Schaefer coordinate CSVs with common column names.

T = readtable(filename, 'VariableNamingRule', 'preserve');
names = T.Properties.VariableNames;

if all(ismember({'R','A','S'}, names))
    CoG = [T.('R'), T.('A'), T.('S')];
elseif all(ismember({'x','y','z'}, names))
    CoG = [T.('x'), T.('y'), T.('z')];
elseif all(ismember({'X','Y','Z'}, names))
    CoG = [T.('X'), T.('Y'), T.('Z')];
else
    error('Unknown coordinate column format in %s.', filename);
end

if ismember('ROI Name', names)
    Labels = T.('ROI Name');
elseif ismember('roi_name', names)
    Labels = T.('roi_name');
else
    Labels = strings(size(CoG,1), 1);
end
end
