function [labels_combined, data_combined, idx_chan_combined] = combine_planar(labels, data, varargin)
% Combines data from planar gradiometers (channels must be first
% dimension). 

parser = inputParser(); 

addParameter(parser, 'method', 'amplitude'); % amplitude or power
addParameter(parser, 'keep_orig_channels', false); 

parse(parser, varargin{:}); 

method = parser.Results.method; 
keep_orig_channels = parser.Results.keep_orig_channels; 

%%

% Magnetometers: Channels ending with 1 (e.g., MEG0111, MEG0121)
% Planar Gradiometers: Channels ending with 2 or 3 (e.g., MEG0112, MEG0113, MEG0122, MEG0123)
grad_cos_all = labels(~cellfun(@isempty, regexp(labels, 'MEG\d+2$')));

assert(size(data, 1) == length(labels), ...
    'Datasize along the first dimension doesnt match the number of channel labels'); 

shape = size(data); 
shape_combined = [length(grad_cos_all), shape(2:end)]; 
data_combined = nan(shape_combined); 
labels_combined = cell(numel(grad_cos_all), 1); 
idx_chan_combined = nan(length(grad_cos_all), 2); 

for i_chan=1:length(grad_cos_all)
    
    % find the orthogonal gradiometer pair
    grad_cos = grad_cos_all{i_chan}; 
    grad_sin = grad_cos; 
    grad_sin(end) = '3'; 
    
    idx_chan = get_chan_idx(labels, {grad_cos, grad_sin}); 
    
    idx_chan_combined(i_chan,:) = idx_chan; 

    assert(length(idx_chan) == 2, ...
        'Failed to find unique grad pair %s + %s', grad_cos, grad_sin)
    
    index = repmat({':'}, 1, ndims(data)); 
    index{1} = idx_chan(1); 
    data_cos = data(index{:}); 
    index{1} = idx_chan(2); 
    data_sin = data(index{:}); 
    
    index{1} = i_chan; 
    if strcmp(method, 'amplitude')
        data_combined(index{:}) = sqrt(data_cos.^2 + data_sin.^2); 
    elseif strcmp(method, 'power')
        data_combined(index{:}) = data_cos + data_sin; 
    else
        error('method ''%s'' not implemented', method)
    end
    
    labels_combined{i_chan} = sprintf('%s+%s', grad_cos, grad_sin(4:end)); 
    
end

if keep_orig_channels
    labels_combined = [labels_combined, labels]; 
    data_combined = cat(1, data_combined, data); 
end

fprintf('%d planars combined assuming data is ''%s'' \n\n', length(labels_combined), method); 

