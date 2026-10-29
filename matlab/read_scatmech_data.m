function [result, source, info] = read_scatmech_data(filename, varargin)
%READ_SCATMECH_DATA Read WaferScope SCATMECH NPZ exports without Python.
% [result,source,info] = read_scatmech_data('far_field_data.npz');
% read_scatmech_data() opens a file chooser. Lists become 1-based cell arrays.
% Numeric axes retain NumPy order; NumPy 1-D vectors become column vectors.
% info.raw is a containers.Map of untouched flat arrays and original keys.
% Requires MATLAB R2016b+ with JVM, no additional toolbox.
if nargin < 1 || isempty(filename)
    [name, folder] = uigetfile('*.npz', 'Select SCATMECH raw data');
    assert(~isequal(name, 0), 'WaferScope:Cancelled', 'No file selected.');
    filename = fullfile(folder, name);
end
parser = inputParser;
addParameter(parser, 'VerifyManifest', true, @(x) islogical(x) && isscalar(x));
parse(parser, varargin{:});
filename = char(filename);
assert(exist(filename, 'file') == 2, 'WaferScope:FileNotFound', ...
    'NPZ file not found: %s', filename);
assert(usejava('jvm'), 'WaferScope:JVMRequired', 'MATLAB JVM is required for ZIP/hash verification.');
fileObject = java.io.File(filename);
filename = char(fileObject.getCanonicalPath());
[folder, stem] = fileparts(filename);
manifestPath = fullfile(folder, [stem '.manifest.json']);
info = struct('filename', filename, 'manifest_verified', [], 'metadata', []);
if parser.Results.VerifyManifest && exist(manifestPath, 'file') == 2
    manifest = read_json(manifestPath);
    actualHash = sha256_file(filename);
    assert(isfield(manifest, 'sha256') && strcmpi(actualHash, manifest.sha256), ...
        'WaferScope:HashMismatch', 'NPZ SHA256 does not match its manifest.');
    info.manifest_verified = true;
    info.sha256 = actualHash;
end
names = checked_zip_entries(filename);
temporary = tempname(tempdir);
mkdir(temporary);
cleaner = onCleanup(@() remove_temporary(temporary)); %#ok<NASGU>
unzip(filename, temporary);
arrays = containers.Map('KeyType', 'char', 'ValueType', 'any');
specifications = containers.Map('KeyType', 'char', 'ValueType', 'any');
for k = 1:numel(names)
    key = names{k}(1:end-4);
    [arrays(key), specifications(key)] = ws_read_npy(fullfile(temporary, names{k}));
end
info.raw = arrays;
info.array_info = specifications;
source = [];
metadataKey = '_archive_metadata_json';
if ~isKey(arrays, metadataKey)
    warning('WaferScope:LegacyArchive', ...
        'Legacy NPZ: flat arrays only; missing nested data cannot be reconstructed.');
    result = struct();
    keys = arrays.keys;
    for k = 1:numel(keys)
        assert(isvarname(keys{k}), 'WaferScope:InvalidKey', ...
            'Legacy key cannot be a MATLAB field; access info.raw instead: %s', keys{k});
        result.(keys{k}) = arrays(keys{k});
    end
    info.mode = 'legacy_flat';
    return
end
text = arrays(metadataKey);
assert(ischar(text) && isrow(text), 'WaferScope:ArchiveMetadata', 'Metadata must be a Unicode scalar.');
validate_json_keys(text);
metadata = jsondecode(text);
assert(strcmp(metadata.format, 'waferscope-raw-result') && isequal(metadata.version, 1), ...
    'WaferScope:ArchiveVersion', 'Unsupported result archive format/version.');
expectedKeys = fieldnames(metadata.arrays);
assert(isequal(sort(arrays.keys), sort([expectedKeys.' {metadataKey}])), ...
    'WaferScope:ArchiveInventory', 'Array inventory does not match metadata.');
for k = 1:numel(expectedKeys)
    key = expectedKeys{k};
    expected = metadata.arrays.(key);
    actual = specifications(key);
    assert(isequal(actual.shape(:), expected.shape(:)) && strcmp(actual.dtype, expected.dtype), ...
        'WaferScope:ArraySpecification', 'Shape or dtype mismatch: %s', key);
end
result = ws_decode_archive(metadata.result, arrays);
if ~isempty(metadata.source)
    source = ws_decode_archive(metadata.source, arrays);
end
info.mode = 'typed_archive_v1';
info.metadata = metadata;
info.metadata_json = text;
end

function names = checked_zip_entries(filename)
archive = java.util.zip.ZipFile(filename);
closer = onCleanup(@() archive.close()); %#ok<NASGU>
entries = archive.entries();
names = {};
while entries.hasMoreElements()
    entry = entries.nextElement();
    name = char(entry.getName());
    assert(~entry.isDirectory() && isempty(regexp(name, '[/\\]', 'once')) && ...
        numel(name) > 4 && strcmp(name(end-3:end), '.npy'), ...
        'WaferScope:UnsafeZipEntry', 'Only flat .npy entries are accepted: %s', name);
    names{end+1} = name; %#ok<AGROW>
end
assert(numel(unique(names)) == numel(names), ...
    'WaferScope:DuplicateArray', 'Duplicate NPZ entries.');
end

function value = read_json(filename)
fid = fopen(filename, 'rb');
assert(fid >= 0, 'WaferScope:JSONOpen', 'Cannot open JSON file.');
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
value = jsondecode(native2unicode(fread(fid, Inf, '*uint8').', 'UTF-8'));
end

function digest = sha256_file(filename)
md = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(filename, 'rb');
assert(fid >= 0, 'WaferScope:HashOpen', 'Cannot open file for SHA256.');
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
while ~feof(fid)
    bytes = fread(fid, 1024^2, '*uint8');
    md.update(typecast(bytes, 'int8'));
end
bytes = typecast(md.digest(), 'uint8');
digest = lower(reshape(dec2hex(bytes, 2).', 1, []));
end

function validate_json_keys(text)
% jsondecode renames invalid identifiers; reject them instead of losing keys.
tokens = regexp(text, '(?<!\\)"(?:[^"\\]|\\.)*"(?=\s*:)', 'match');
for k = 1:numel(tokens)
    key = jsondecode(tokens{k});
    assert(isvarname(key), 'WaferScope:InvalidKey', ...
        'JSON key is not a MATLAB identifier; refusing to rename it: %s', key);
end
end

function remove_temporary(folder)
% Delete only the unique extraction directory created under tempdir.
parentObject = java.io.File(fileparts(folder));
tempObject = java.io.File(tempdir);
parent = char(parentObject.getCanonicalPath());
expected = char(tempObject.getCanonicalPath());
if strcmp(parent, expected) && exist(folder, 'dir') == 7
    rmdir(folder, 's');
end
end
