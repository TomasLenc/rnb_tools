function [labels_combined, idx_chan_combined] = combine_planar_labels(labels)

grad_cos_all = labels(~cellfun(@isempty, regexp(labels, 'MEG\d+2$')));

labels_combined = cell(size(grad_cos_all)); 
idx_chan_combined = nan(length(grad_cos_all), 2); 

for i_chan=1:length(grad_cos_all)
    
    % find the orthogonal gradiometer pair
    grad_cos = grad_cos_all{i_chan}; 
    grad_sin = grad_cos; 
    grad_sin(end) = '3'; 
    
    idx_chan = get_chan_idx(labels, {grad_cos, grad_sin}); 
    
    assert(length(idx_chan) == 2, ...
        'Failed to find unique grad pair %s + %s', grad_cos, grad_sin)
    
    labels_combined{i_chan} = sprintf('%s+%s', grad_cos, grad_sin(4:end)); 
    idx_chan_combined(i_chan,:) = idx_chan; 
end