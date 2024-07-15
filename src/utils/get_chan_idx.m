function idx = get_chan_idx(header_chanlocs_or_labels, chan_labels)

% determine if this is header, chanlocs, or labels
if iscell(header_chanlocs_or_labels)
    % we got a cell of label strings
    labs = header_chanlocs_or_labels; 
elseif isstruct(header_chanlocs_or_labels) && isfield(header_chanlocs_or_labels, 'chanlocs')
    % we got a header 
    labs = {header_chanlocs_or_labels.chanlocs.labels}; 
elseif isstruct(header_chanlocs_or_labels) && isfield(header_chanlocs_or_labels, 'labels')
    % we got a chanlocs structure
    labs = {header_chanlocs_or_labels.labels}; 
else
    error('cannot figure out what the input is...'); 
end

idx = cellfun(@(x)find(strcmpi(labs, x)),  chan_labels, 'uni', 1);

