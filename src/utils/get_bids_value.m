function value = get_bids_value(fname, key)

% tmp = regexp(fname, sprintf('%s-([^_]+)(?:_|$)', key), 'tokens', 'once');
tmp = regexp(fname, sprintf('%s-([^_./]+)', key), 'tokens', 'once'); 

if ~isempty(tmp)
    value = tmp{1};
else
    value = [];
end