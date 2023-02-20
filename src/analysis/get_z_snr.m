function [z_snr, mean_snip, idx_snip] = get_z_snr(mX, freq, frex, bin_from, bin_to)
% Calculate SNR as zscore across harmonics. 
% 
% Parameters
% ----------
% mX : array_like, shape=[..., frequency]
%     Raw (not noise-subtracted!!!) magnitude spectra with frequency as the 
%     last dimension. 
% freq : array_like
%     Frequencies for the FFT. 
% frex : array_like
%     Frequencies to take the snippets around. 
% bin_from : int
%     Closest frequency bin (on both sides) used to calculate noise properties.
% bin_to : int
%     Furthest frequency bin (on both sides) used to calculate noise properties.
% 
% Returns 
% -------
% z_snr : array_like, shape=shape(mX)(1:end-1)
%     Z-score values from the center bin taken from average snippet around all 
%     frex. 
% mean_snip : array_like, shape=shape(mX)(1:end-1)
%     Average snippet around all frex. 
% idx_snip : array_like
%     Relative index of each bin in the mean_snip

relative_idx = [-bin_to : +bin_to]; 

shape = size(mX); 
snippets = nan([shape(1:end-1), length(relative_idx), length(frex)]); 

N = size(mX, ndims(mX)); 

frex_idx = dsearchn(freq', frex'); 

for i_f=1:length(frex)
    
    idx_start = frex_idx(i_f) - bin_to; 
    idx_end = frex_idx(i_f) + bin_to; 
    
    if idx_start <= 0
        error('Frex %d is too close to the start of the input (idx_start = %d).', ...
              i_f, idx_start); 
    end
    if idx_end > shape(end)
        error('Frex %d is too close to the end of the input (idx_end = %d).', ...
              i_f, idx_end); 
    end

    index_snip = repmat({':'}, 1, ndims(mX) + 1);
    index_snip{ndims(mX) + 1} = i_f;
    
    index_mX = repmat({':'}, 1, ndims(mX));
    index_mX{end} = [idx_start : idx_end];

    snippets(index_snip{:}) = mX(index_mX{:}); 

end

mean_snip = mean(snippets, ndims(mX) + 1); 

idx_noise = [-bin_to : -bin_from, +bin_from : +bin_to]; 
idx_noise = idx_noise - idx_noise(1) + 1; 

idx_signal = bin_to + 1; 

index = repmat({':'}, 1, ndims(mX));
index{end} = idx_signal; 
s_mean = mean_snip(index{:}); 

index = repmat({':'}, 1, ndims(mX));
index{end} = idx_noise; 
noise_mean = mean(mean_snip(index{:}), ndims(mX)); 

index = repmat({':'}, 1, ndims(mX));
index{end} = idx_noise; 
noise_sd = std(mean_snip(index{:}), [], ndims(mX)); 

z_snr = (s_mean - noise_mean) ./ noise_sd; 

% figure
% plot(mean_snip)
% hold on
% plot(idx_noise, mean_snip(idx_noise), 'bo')
% plot(idx_signal, mean_snip(idx_signal), 'ro')

% this is just for output
idx_snip = [-bin_to : +bin_to]; 