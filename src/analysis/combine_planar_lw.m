function [header_combined, data_combined] = combine_planar_lw(header, data, varargin)
% Combines planar gradiometers. Input is lw6 data format. 

parser = inputParser(); 

addParameter(parser, 'method', 'amplitude'); % amplitude or power
addParameter(parser, 'keep_orig_channels', true); 

parse(parser, varargin{:}); 

method = parser.Results.method; 
keep_orig_channels = parser.Results.keep_orig_channels; 

%%

% Magnetometers: Channels ending with 1 (e.g., MEG0111, MEG0121)
% Planar Gradiometers: Channels ending with 2 or 3 (e.g., MEG0112, MEG0113, MEG0122, MEG0123)


channels = {header.chanlocs.labels}; 

grad_cos_all = channels(~cellfun(@isempty, regexp(channels, 'MEG\d+2$')));

data_combined = nan(header.datasize(1), length(grad_cos_all), header.datasize(3), header.datasize(4), header.datasize(5), header.datasize(6)); 
chanlocs_combined = []; 

for i_chan=1:length(grad_cos_all)
    
    % find the orthogonal gradiometer pair
    grad_cos = grad_cos_all{i_chan}; 
    grad_sin = grad_cos; 
    grad_sin(end) = '3'; 
    
    idx_chan = get_chan_idx(header, {grad_cos, grad_sin}); 
    
    assert(length(idx_chan) == 2, ...
        'Failed to find unique grad pair %s + %s', grad_cos, grad_sin)
    
    if strcmp(method, 'amplitude')
        data_combined(:, i_chan, :, :, :, :) = sqrt(data(:, idx_chan(1), :, :, :, :).^2 + data(:, idx_chan(2), :, :, :, :).^2); 
    elseif strcmp(method, 'power')
        data_combined(:, i_chan, :, :, :, :) = data(:, idx_chan(1), :, :, :, :) + data(:, idx_chan(2), :, :, :, :); 
    else
        error('method ''%s'' not implemented', method)
    end
    
    chanlocs_combined(i_chan).labels = sprintf('%s+%s', grad_cos, grad_sin(4:end)); 
    chanlocs_combined(i_chan).topo_enabled = false; 
    chanlocs_combined(i_chan).SEEG_enabled = false; 
    
end

if keep_orig_channels
    chanlocs_combined = [chanlocs_combined, header.chanlocs]; 
    data_combined = cat(2, data_combined, data); 
else
    mask_grad = ~cellfun(@isempty, regexp(channels, 'MEG\d+[2,3]$'));
    chanlocs_combined = [chanlocs_combined, header.chanlocs(~mask_grad)]; 
    data_combined = cat(2, data_combined, data(:, ~mask_grad, :, :, :, :)); 
end

header_combined = header; 
header_combined.datasize = size(data_combined); 
header_combined.chanlocs = chanlocs_combined; 

fprintf('%d planars combined assuming data is ''%s'' \n\n', length(grad_cos_all), method); 

