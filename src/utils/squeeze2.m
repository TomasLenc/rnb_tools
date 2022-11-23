function out = squeeze2(A, dims_to_squeeze)
% This is an extension of the default squeeze() function, where you can
% specify which dimensions are going to be squeezed. The rest will be
% untouched. 

sz = size(A); 

% make sure that all out_dim_mask are actually singleton in the input
% array
assert(all(sz(dims_to_squeeze) == 1))

% make a binary mask about which dimensions are non-singleton
out_dim_mask = true(1, ndims(A)); 

% set out_dim_mask to false in the mask  
out_dim_mask(dims_to_squeeze) = false; 

out = reshape(A, sz(out_dim_mask)); 
