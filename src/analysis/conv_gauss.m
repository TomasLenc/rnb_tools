function [x_out, fs_decim] = conv_gauss(x, fs, varargin)
% Convolve time-domain signal x with a gaussian kernel with given
% parameters. 
%
% Parameters
% ----------
% x : shape=[1, time]
%     Input time series. 
% fs : int
%     Sampling rate.
% decim_factor : int, optional, default=1
%     Optional decimation of the signal after convolution. 
% fwhm_time : float, optional, default=[]
%     Full Width at Half Maximum (time-domain) of the gaussian kernel in
%     seconds.
% sd_time : float, optional, default=[]
%     Standard deviation of the gaussian kernel in seconds. 
% 
% Returns 
% -------
% x_out : shape=[1, time]
%     Output time series. 
% fs_decim : float
%     Sampling rate of the output signal. 

parser = inputParser; 

addParameter(parser, 'decim_factor', 1); 
addParameter(parser, 'fwhm_time', []); 
addParameter(parser, 'sd_time', []); 

parse(parser, varargin{:}); 

decim_factor = parser.Results.decim_factor; 
fwhm = parser.Results.fwhm_time; 
sigma = parser.Results.sd_time; 

assert(~isempty(fwhm) || ~isempty(sigma), ...
       'Please provide either FWHM or SD value for the gaussian kernel!')

%%

if ~isempty(fwhm) 
    
    t = [-(fwhm * 10) : 1/fs : +(fwhm * 10)]; 
    
    kernel = exp( (-4 * log(2) * t.^2) ./ fwhm^2 );

elseif ~isempty(sigma) 
    
    t = [-(sigma * 10) : 1/fs : +(sigma * 10)]; 
    
    kernel = 1 / sqrt(2*pi*sigma^2) * exp(-1/2 * (t ./ sigma).^2);
    % y = normpdf(t, 0, sigma); 

end

[x_out, fs_decim] = convolve_fft(x, kernel, fs,...
                                 'decim_factor', decim_factor); 

