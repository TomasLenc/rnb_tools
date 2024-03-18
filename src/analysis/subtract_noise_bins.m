function mX_subtracted = subtract_noise_bins(mX, from, to)
% Go over the magnitude spectrum and for each frequency bin, subtract the
% mean taken across surrounding bins. 
% 
% Parameters
% ----------
% mX : array_like, shape=[..., frequency]
%     Raw magnitude spectra with frequency as the last dimension.
% from : int
%     Closest frequency bin (on both sides) used to calculate noise properties.
% to : int
%     Furthest frequency bin (on both sides) used to calculate noise properties.
% 
% Returns 
% -------
% mX_subtracted : array_like, shape=shape(mX)
%     Noise-subtracted magnitude spectrum. 
%

mX_subtracted = mX; 

N = size(mX, ndims(mX)); 

for i_f=1:length(mX)
    
    idx_1 = max(i_f - to, 1); 
    idx_2 = max(i_f - from, 1); 
    idx_3 = min(i_f + from, N); 
    idx_4 = min(i_f + to, N); 

    index = cell(1, ndims(mX));
    index(:) = {':'};
    index{end} = [idx_1:idx_2, idx_3:idx_4];
    mean_noise = mean(mX(index{:}), ndims(mX)); 

    index = cell(1, ndims(mX));
    index(:) = {':'};
    index{end} = i_f;
    mX_subtracted(index{:}) = mX(index{:}) - mean_noise; 
end