function X_subtracted = subtract_noise_bins_complex(X, from, to)
% Go over the complex spectrum and for each frequency bin, prepare a noise
% vector that has the same direction to the complex vector [real, imag] at
% the center bin, and its magnitude is equal to the average magnitude of
% complex vectors taken from surrounding bins. Then subtract this noise
% vector from the vector at the center bin. If the noise vector has larger
% magnitude than the vector at the center bin, set the resulting complex
% vector to be the zero vector. This prevent overshooting and getting a big
% vector with opposite direction to the original vector at the center bin.
% 
% Parameters
% ----------
% X : array_like, shape=[..., frequency]
%     Raw complex DFT spectra with frequency as the last dimension.
% from : int
%     Closest frequency bin (on both sides) used to calculate noise properties.
% to : int
%     Furthest frequency bin (on both sides) used to calculate noise properties.
% 
% Returns 
% -------
% X_subtracted : array_like, shape=shape(mX)
%     Noise-subtracted complex DFT spectrum. 
%

mX = abs(X); 

X_subtracted = X; 

N = size(X, ndims(X)); 

for i_f=1:length(mX)
    
    idx_1 = max(i_f - to, 1); 
    idx_2 = max(i_f - from, 1); 
    idx_3 = min(i_f + from, N); 
    idx_4 = min(i_f + to, N); 

    index = cell(1, ndims(mX));
    index(:) = {':'};
    index{end} = [idx_1:idx_2, idx_3:idx_4];
    
    mean_noise_magn = mean(mX(index{:}), ndims(mX)); 

    index = cell(1, ndims(mX));
    index(:) = {':'};
    index{end} = i_f;
    
    % get the vector that we will add to the vector at the target bin
    mean_noise_vect = X(index{:}) ./ abs(X(index{:})); 
    mean_noise_vect = - mean_noise_vect .* mean_noise_magn; 
    
%     % sanity check figure
%     figure
%     plot([0, real(mean_noise_vect(1))], [0, imag(mean_noise_vect(1))], '-', 'linew', 3)
%     hold on 
%     plot([0, real(X(1, i_f))], [0, imag(X(1, i_f))], '-', 'linew', 3)
%     ax = gca; 
%     ax.XAxisLocation = 'origin'; 
%     ax.YAxisLocation = 'origin'; 
    
    vals_subtr = X(index{:}) + mean_noise_vect; 
    
    % check where the noise magnitude is larger than the magnitude at the
    % target bin and set the resulting vector at the target bin to zero
    mask_to_zero = abs(mean_noise_vect) > abs(X(index{:})); 
    vals_subtr(mask_to_zero) = 0; 
    
    X_subtracted(index{:}) = vals_subtr; 
    
end