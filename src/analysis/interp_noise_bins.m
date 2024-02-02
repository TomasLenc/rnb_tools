function mX_interp = interp_noise_bins(mX, idx_to_interp, from, to)
% Replace magnitudes at the requested frequency bins with mean of the bins
% around. Frequecy must be the last dimension. 

mX_interp = mX; 

N = size(mX, ndims(mX)); 

for i_f=1:length(idx_to_interp)
    idx_1 = max(idx_to_interp(i_f) - to, 1); 
    idx_2 = max(idx_to_interp(i_f) - from, 1); 
    idx_3 = min(idx_to_interp(i_f) + from, N); 
    idx_4 = min(idx_to_interp(i_f) + to, N); 

    index = cell(1, ndims(mX));
    index(:) = {':'};
    index{end} = [idx_1:idx_2, idx_3:idx_4];
    mean_around_bin = mean(mX(index{:}), ndims(mX)); 

    index = cell(1, ndims(mX));
    index(:) = {':'};
    index{end} = idx_to_interp(i_f);
    mX_interp(index{:}) = mean_around_bin; 
end
