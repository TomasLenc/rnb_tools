function idx = get_chan_idx(header_chanlocs_or_labels, chan_labels, varargin)

parser = inputParser;

addParameter(parser, 'nan_not_found', true);

parse(parser, varargin{:});

nan_not_found = parser.Results.nan_not_found;

%%

% convert to cell if necessary
if ischar(chan_labels)
    chan_labels = {chan_labels}; 
end

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

% go over requested channel labels and find each 
idx = nan(1, length(chan_labels)); 

for i_chan=1:length(chan_labels)
    
    tmp = find(strcmpi(labs, chan_labels{i_chan})); 
    
    if isempty(tmp)
        warning('chan %s found no matches', chan_labels{i_chan}); 
    elseif length(tmp) > 1
        warning('chan %s found more than 1 match', chan_labels{i_chan}); 
    else
        idx(i_chan) = tmp; 
    end
    
end

if ~nan_not_found
    idx = idx(~isnan(idx));
end

