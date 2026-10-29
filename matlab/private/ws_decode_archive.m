function value = ws_decode_archive(node, arrays)
%WS_DECODE_ARCHIVE Resolve typed nodes by their saved array keys, not names.
assert(isstruct(node) && isscalar(node) && isfield(node, 'kind'), ...
    'WaferScope:ArchiveNode', 'Invalid archive tree node.');
switch node.kind
    case 'array'
        assert(isKey(arrays, node.key), 'WaferScope:MissingArray', ...
            'Missing array referenced by metadata: %s', node.key);
        value = arrays(node.key);
    case 'dict'
        value = struct();
        names = fieldnames(node.items);
        for k = 1:numel(names)
            value.(names{k}) = ws_decode_archive(node.items.(names{k}), arrays);
        end
    case {'list', 'tuple'}
        items = node.items;
        value = cell(1, numel(items));
        for k = 1:numel(items)
            if iscell(items)
                item = items{k};
            else
                item = items(k);
            end
            value{k} = ws_decode_archive(item, arrays);
        end
    case 'complex'
        value = complex(node.real, node.imag);
    case 'scalar'
        value = node.value;
    case 'nonfinite'
        switch node.value
            case 'nan'
                value = NaN;
            case 'inf'
                value = Inf;
            case '-inf'
                value = -Inf;
            otherwise
                error('WaferScope:ArchiveNode', 'Unknown nonfinite value.');
        end
    otherwise
        error('WaferScope:ArchiveNode', 'Unknown node kind: %s.', node.kind);
end
end
