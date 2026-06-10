function HE_out = unique_hyperedges(HE_in)

if isempty(HE_in)
    HE_out = {};
    return;
end

tmp = cell(size(HE_in));
for i = 1:numel(HE_in)
    tmp{i} = sort(unique(HE_in{i}));
end

keys = cellfun(@(x) sprintf('%d_', x), tmp, 'UniformOutput', false);
[~, ia] = unique(keys, 'stable');
HE_out = tmp(ia);

end