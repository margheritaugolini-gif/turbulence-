function data = readNPY(filename)
%READNPY Read a simple numeric NumPy .npy file into MATLAB.
% Supports little-endian numeric arrays stored in .npy v1.0/v2.0 format.

fid = fopen(filename, 'rb');
if fid < 0
    error('Cannot open NPY file: %s', filename);
end
cleanup = onCleanup(@() fclose(fid));

magic = fread(fid, 6, '*uint8')';
if ~isequal(magic, uint8([147 double('NUMPY')]))
    error('File is not a NumPy .npy file: %s', filename);
end

version = fread(fid, 2, '*uint8')';
if version(1) == 1
    headerLen = fread(fid, 1, 'uint16', 0, 'ieee-le');
elseif version(1) == 2
    headerLen = fread(fid, 1, 'uint32', 0, 'ieee-le');
else
    error('Unsupported NPY version %d.%d.', version(1), version(2));
end

header = char(fread(fid, headerLen, '*char')');

descr = regexp(header, '''descr'':\s*''([^'']+)''', 'tokens', 'once');
fortranOrder = regexp(header, '''fortran_order'':\s*(True|False)', 'tokens', 'once');
shapeText = regexp(header, '''shape'':\s*\(([^)]*)\)', 'tokens', 'once');

if isempty(descr) || isempty(fortranOrder) || isempty(shapeText)
    error('Could not parse NPY header in %s.', filename);
end

[precision, machinefmt] = npyPrecision(descr{1});
shape = sscanf(shapeText{1}, '%d,')';
if isempty(shape)
    error('Could not parse NPY shape in %s.', filename);
end

count = prod(shape);
raw = fread(fid, count, ['*' precision], 0, machinefmt);
if numel(raw) ~= count
    error('Expected %d values in %s, found %d.', count, filename, numel(raw));
end

if strcmp(fortranOrder{1}, 'True')
    data = reshape(raw, shape);
elseif numel(shape) == 2
    data = reshape(raw, fliplr(shape))';
else
    error('C-order NPY arrays with more than 2 dimensions are not supported.');
end
end

function [precision, machinefmt] = npyPrecision(descr)
endian = descr(1);
dtype = descr(2:end);

switch endian
    case '<'
        machinefmt = 'ieee-le';
    case '>'
        machinefmt = 'ieee-be';
    case '|'
        machinefmt = 'native';
    otherwise
        error('Unsupported NPY endian marker: %s', endian);
end

switch dtype
    case 'f4'
        precision = 'single';
    case 'f8'
        precision = 'double';
    case 'i4'
        precision = 'int32';
    case 'i8'
        precision = 'int64';
    case 'u4'
        precision = 'uint32';
    case 'u8'
        precision = 'uint64';
    case 'i2'
        precision = 'int16';
    case 'u2'
        precision = 'uint16';
    case 'i1'
        precision = 'int8';
    case 'u1'
        precision = 'uint8';
    otherwise
        error('Unsupported NPY dtype: %s', descr);
end
end
