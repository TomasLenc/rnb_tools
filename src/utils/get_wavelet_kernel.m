function [cmw, fig_ir] = get_wavelet_kernel(f, sigma, fs, varargin)
% Prepare complex morlet wavelet in the time domain, based on the requested
% center frequency and standard deviation of the gaussian in the frequency
% domain. 
%
% Parameters
% ----------
% f : float
%     Center frequency in Hz. 
% sigma : float
%     Standard deviation in the frequency domain (Hz). 
% fs : float
%     Sampling rate (samples/s); 
% do_plot: bool, optional, default=false
%     If true, diagnostic plots of the kernel are generated, with empirical 
%     FWHM. 
%     
% Returns
% -------
% cmw : array of complex floats
%     Complex morlet wavelet kernel. 

% f = 1; 
% sigma = f/2; 
% fs = 44100; 

parser = inputParser; 

addParameter(parser, 'do_plot', false); 

parse(parser, varargin{:});

do_plot = parser.Results.do_plot; 


%%

t = [-10: 1/fs: 10]; 
n = length(t); 
n_conv_pow2 = 2 ^ nextpow2(round(fs * 10)); 
freq = [0 : floor(n_conv_pow2/2)] / n_conv_pow2 * fs; 

% we need to make the spectra column-wise because of bsxfun later...trasposing this
% later would do complex conjugate!
    
% parametrization based on frequency-domain SD value
A = 1 / sqrt(sigma * sqrt(pi)); 
gwin = A * exp(- t.^2 / (2 * sigma^2)); 

% % parametrization based on time-domain FWHM value
% gwin = exp( (-4*log(2)*wt.^2) ./ fwhmT^2 );

% % Mike Cohen has a mistake in his code, squaring the whole thing...
% gwin = exp(-(4*log(2)*wt).^2/fwhmT.^2)
cos_wave = exp( 1j * 2*pi * f * t ); 

cmw = cos_wave .* gwin;    
    
%% empirical wavelet properties

% magnitude spectrum of the wavelet
cmwX = fft(cmw, n_conv_pow2); 

% Normalization of the wavelet for accurate reconstruction of the signal 
% amplitude can be achieved by setting the amplitude of the Fourier transform 
% of the wavelet to have a peak of 1; this is analogous to conceptualizing the 
% frequency response of a temporal filter as a gain function and setting the 
% passband frequencies with a gain of 1. 
cmwX = cmwX ./ max(cmwX);

% compute the empirical temporal FWHM in seconds (later converted to ms)
gwin_norm = gwin / max(gwin); % normalize to 1
midp = dsearchn(ensure_col(t), 0);
fwhmT_upper = t(midp - 1 + dsearchn(ensure_col(gwin_norm(midp:end)), .5)); 
fwhmT_lower = t(dsearchn(ensure_col(gwin_norm(1:midp)), .5)); 
fwhmT = fwhmT_upper - fwhmT_lower;

% compute the empirical frequency FWHM in Hz
cmwmX_norm  = cmwX / max(cmwX); % normalize to 1 to get the right units 
% cmwX_norm = cmwX_norm.^2; % get power (not amplitude)??? (I saw both in Cohen...sometimes he just calls amplitude power...or maybe he has a mistake in the code?)
frex_idx = dsearchn(ensure_col(freq), f);
fwhmF_upper = freq(frex_idx - 1 + dsearchn(ensure_col(abs(cmwmX_norm(frex_idx:end))), .5)); 
fwhmF_lower = freq(dsearchn(ensure_col(abs(cmwmX_norm(1:frex_idx))), .5)); 
fwhmF = fwhmF_upper - fwhmF_lower;

%% plot 

if do_plot
    
    fig_ir = figure;

    ax = subplot(2, 1, 1); 
    cmw_norm = cmw / max(gwin); 
    plot(t(1:1:end), real(cmw_norm(1:1:end))); 
    hold on
    plot([fwhmT_lower, fwhmT_lower], [-1,1], 'r:'); 
    plot([fwhmT_upper, fwhmT_upper], [-1,1], 'r:'); 
    title(sprintf('f=%.1fHz \nfwhmT=%.3fs', f, fwhmT)); 

    ax = subplot(2, 1, 2); 
    cmw_plot_maxfreq = min(fs/2, f*2); 
    cmw_plot_maxfreq_idx = round(cmw_plot_maxfreq / fs * n_conv_pow2) + 1; 
    cmw_freq_plot = [0 : cmw_plot_maxfreq_idx-1] / n_conv_pow2 * fs; 
    plot(cmw_freq_plot, abs(cmwmX_norm(1:length(cmw_freq_plot)))); 
    hold on
    plot([fwhmF_lower, fwhmF_lower], [0,1], 'r:'); 
    plot([fwhmF_upper, fwhmF_upper], [0,1], 'r:'); 
    title(sprintf('fwhmF=%.2fHz', fwhmF)); 
    ax.XLim = [0, cmw_plot_maxfreq]; 
    ax.XTick = [f - f/2, f, f + f/2]; 
    ax.YTick = []; 

    fprintf('f = %.1f Hz \t fwhmT = %.3f sec \t fwhmF = %.3f Hz \n', ...
            f, fwhmT, fwhmF); 

else
   fig_ir = [];  
end


