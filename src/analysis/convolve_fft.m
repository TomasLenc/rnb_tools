function [x_out, fs_decim] = convolve_fft(x, kernel, fs, varargin)
% Convolution of a signal x with a kernel in the frequency domain (based on
% Cohen 2014). Optionally also downsamples the data (using matlab function
% decimate, which applies low-pass filter before downsampling). 

parser = inputParser; 

addParameter(parser, 'decim_factor', 1); 

parse(parser, varargin{:}); 

decim_factor = parser.Results.decim_factor; 

%%

fs_decim = fs / decim_factor; 

n_data          = length(x); 
x_dur           = length(x) / fs; 
n_out_data      = round(x_dur * fs_decim); 

n_w             = length(kernel); 
n_hw            = floor(n_w / 2); 
n_conv          = n_data + n_w - 1; 
n_conv_pow2     = pow2(nextpow2(n_conv)); 

%%

% FFT of data
X = fft(x, n_conv_pow2, 2); 

% FFT of kernel
X_kernel = fft(kernel, n_conv_pow2); 

% convolution 
as = abs(ifft(bsxfun(@times, X_kernel.', X.'), [], 1)).'; 

% cut-off edges
as = as(:, n_hw+1 : n_conv-n_hw); 

%%

% downsample
if decim_factor ~= 1
    x_out = nan(size(as, 1), n_out_data); 
    for i_trial=1:size(as, 1)
        tmp = decimate(as(i_trial, :), decim_factor); 
        x_out = tmp(1:n_out_data); 
    end
else
    x_out = as; 
end
