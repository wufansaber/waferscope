function [value, info] = ws_read_npy(filename)
%WS_READ_NPY Read primitive NPY arrays without Python or pickle.
fid = fopen(filename, 'rb');
assert(fid >= 0, 'WaferScope:NPYOpen', 'Cannot open NPY: %s', filename);
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
magic = fread(fid, 6, '*uint8').';
assert(isequal(magic, uint8([147 double('NUMPY')])), ...
    'WaferScope:NPYMagic', 'Invalid NPY signature.');
version = fread(fid, 2, '*uint8');
assert(numel(version) == 2 && version(2) == 0 && any(version(1) == [1 2 3]), ...
    'WaferScope:NPYVersion', 'Unsupported NPY version.');
if version(1) == 1
    headerLength = fread(fid, 1, 'uint16=>double', 0, 'ieee-le');
else
    headerLength = fread(fid, 1, 'uint32=>double', 0, 'ieee-le');
end
assert(~isempty(headerLength) && headerLength <= 1024^2, ...
    'WaferScope:NPYHeader', 'Invalid or oversized NPY header.');
header = char(read_exact(fid, headerLength, 'uint8', 'native').');
dtype = regexp(header, '''descr''\s*:\s*''([^'']+)''', 'tokens', 'once');
order = regexp(header, '''fortran_order''\s*:\s*(True|False)', 'tokens', 'once');
shape = regexp(header, '''shape''\s*:\s*\(([^)]*)\)', 'tokens', 'once');
assert(~isempty(dtype) && ~isempty(order) && ~isempty(shape), ...
    'WaferScope:NPYHeader', 'Unsupported NPY header or structured dtype.');
shapeText = strtrim(shape{1});
if isempty(shapeText)
    dimensions = [];
else
    assert(~isempty(regexp(shapeText, '^(?:[0-9]+\s*,\s*)*(?:[0-9]+\s*)?$', 'once')), ...
        'WaferScope:NPYShape', 'Invalid NPY dimensions.');
    dimensions = str2double(regexp(shapeText, '[0-9]+', 'match'));
end
count = prod(dimensions);
assert(all(isfinite(dimensions)) && all(dimensions <= flintmax) && count <= flintmax, ...
    'WaferScope:NPYShape', 'NPY dimensions exceed exact indexing range.');
type = regexp(dtype{1}, '^([<>=|])([fiucbU])([0-9]+)$', 'tokens', 'once');
assert(~isempty(type), 'WaferScope:NPYDtype', ...
    'Unsupported dtype %s; object/pickle arrays are never loaded.', dtype{1});
width = str2double(type{3});
machine = 'native';
if strcmp(type{1}, '<')
    machine = 'ieee-le';
elseif strcmp(type{1}, '>')
    machine = 'ieee-be';
end
itemBytes = width;
if strcmp(type{2}, 'U')
    itemBytes = 4*width;
end
position = ftell(fid);
fseek(fid, 0, 'eof');
assert(ftell(fid)-position == count*itemBytes, ...
    'WaferScope:NPYLength', 'NPY payload size does not match shape/dtype.');
fseek(fid, position, 'bof');
if strcmp(type{2}, 'U')
    codes = read_exact(fid, count*width, 'uint32', machine);
    codes = reshape(codes, width, count);
    text = cell(count, 1);
    for k = 1:count
        text{k} = utf32_text(codes(:, k));
    end
    if isempty(dimensions)
        value = text{1};
    else
        value = numpy_shape(text, dimensions, strcmp(order{1}, 'True'));
    end
else
    precision = numeric_precision(type{2}, width);
    if strcmp(type{2}, 'c')
        parts = read_exact(fid, 2*count, precision, machine);
        flat = complex(parts(1:2:end), parts(2:2:end));
    else
        flat = read_exact(fid, count, precision, machine);
        if strcmp(type{2}, 'b')
            flat = logical(flat);
        end
    end
    value = numpy_shape(flat, dimensions, strcmp(order{1}, 'True'));
end
info = struct('shape', dimensions, 'dtype', dtype{1}, ...
    'fortran_order', strcmp(order{1}, 'True'));
end

function value = read_exact(fid, count, precision, machine)
value = fread(fid, count, ['*' precision], 0, machine);
assert(numel(value) == count, 'WaferScope:NPYTruncated', 'Truncated NPY data.');
end

function precision = numeric_precision(kind, width)
if strcmp(kind, 'f') && any(width == [4 8])
    names = {'single', 'double'};
    precision = names{1+(width == 8)};
elseif strcmp(kind, 'c') && any(width == [8 16])
    names = {'single', 'double'};
    precision = names{1+(width == 16)};
elseif any(strcmp(kind, {'i', 'u'})) && any(width == [1 2 4 8])
    prefix = 'int';
    if strcmp(kind, 'u')
        prefix = 'uint';
    end
    precision = [prefix num2str(8*width)];
elseif strcmp(kind, 'b') && width == 1
    precision = 'uint8';
else
    error('WaferScope:NPYDtype', 'Unsupported numeric NPY dtype: %s%d.', kind, width);
end
end

function value = numpy_shape(flat, dimensions, fortranOrder)
if isempty(dimensions)
    value = flat(1);
elseif numel(dimensions) == 1
    value = reshape(flat, dimensions(1), 1);
elseif fortranOrder
    value = reshape(flat, dimensions);
else
    value = permute(reshape(flat, fliplr(dimensions)), numel(dimensions):-1:1);
end
end

function text = utf32_text(codes)
codes = double(codes(:).');
last = find(codes ~= 0, 1, 'last');
if isempty(last)
    text = '';
    return
end
codes = codes(1:last);
assert(all(codes <= hex2dec('10FFFF')) && ...
    ~any(codes >= hex2dec('D800') & codes <= hex2dec('DFFF')), ...
    'WaferScope:NPYUnicode', 'Invalid UTF-32 code point.');
wide = codes > 65535;
lengths = 1+wide;
positions = [1 1+cumsum(lengths(1:end-1))];
units = zeros(1, sum(lengths), 'uint16');
units(positions(~wide)) = uint16(codes(~wide));
extra = codes(wide)-65536;
units(positions(wide)) = uint16(55296+floor(extra/1024));
units(positions(wide)+1) = uint16(56320+mod(extra, 1024));
text = char(units);
end
