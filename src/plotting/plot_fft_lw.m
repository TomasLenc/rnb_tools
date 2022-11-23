function plot_fft_lw(header, data, varargin)
% This assumes .lw6 files. No FFT or SNR subtraction is performed here, we
% assume data already containts spectra that are ready for plotting. 

freq = [0 : header.datasize(end) - 1] * header.xstep; 

mX = squeeze(data)'; 
mX(1) = 0; 

plot_fft(freq, mX, varargin{:}); 