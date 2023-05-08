function x = add_signal_noise(x_clean, noise, snr)
% Add signal and noise timeseries to yield a requested SNR. 
% 
% Parameters
% ----------
% x_clean : array_like
%     Clean signal timeseries. Time must be the last dimension. 
% x_clean : array_like
%     Noise timeseries to be mixed with the signal. Must have the same size 
%     as x_clean. 
% snr : float
%     Requested signal-to-noise ratio. 
% 
% Returns
% -------
% x : array_like
%     Mixed signal and noise with the requested SNR. 

x_clean_rms = rms(x_clean); 
noise_rms = rms(noise, ndims(noise)); 
noise_gain = (x_clean_rms ./ noise_rms) / snr; 
noise = noise .* noise_gain; 

% add signal and noise
x = x_clean + noise; 