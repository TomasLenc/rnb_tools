function x = get_colored_noise(N, fs, exponent, varargin)
% Genrates 1-D colored noise. 
% 
% Parameters
% ----------
% N : int
%     Number of requested samples. 
% fs : int
%     Sampling rate.
% exponent : float
%     Exponent f^beta (careful! not 1/(f^beta)) so don't mess up the sign ;)
% method : {'cohen', 'neurodsp'}
%     Which method to use for synthesis. 
% f_rotation : float
%     If `method` is 'neurodsp', this is the frequency around which the spectrum
%     will be rotated. 
% 
% Returns
% -------
% x : array_like
%     Colored noise with the requested exponent. 

parser = inputParser; 

addParameter(parser, 'f_rotation', 1); 
addParameter(parser, 'method', 'cohen'); 

parse(parser, varargin{:});
f_rotation = parser.Results.f_rotation; 
method = parser.Results.method; 


if strcmp(method, 'cohen')

    % method from Cohen 

    T = N/fs; 
    freq = [0:N-1]/T; 

    if mod(N,2)==0
        hN          = N/2-1; 
        magnitudes  = rand(1,hN) .* (freq(2:hN+1) .^ (exponent/2)); 
        phases      = exp(1j*2*pi*rand(1,hN)); 
        X           = [0, magnitudes, 0, magnitudes(end:-1:1)] .* [0, phases, 0, conj(phases(end:-1:1))]; 
    else
        hN          = floor(N/2); 
        magnitudes  = rand(1,hN) .* (freq(2:hN+1) .^ (exponent/2)); 
        phases      = exp(1j*2*pi*rand(1,hN)); 
        X           = [0, magnitudes, magnitudes(end:-1:1)] .* [0, phases, conj(phases(end:-1:1))]; 
    end

    x = real(ifft(X)); 
    
elseif strcmp(method, 'neurodsp')

    % method from neurodsp
    x = randn(1, N); 
    
    T = N/fs; 
    freq = [0:N-1]/T; 
    
    % Delta exponent is divided by two, as the FFT output is in units of 
    % amplitude not power. 
    delta_exponent = exponent / 2; 
    
    X = fft(x); 
    
    if freq(1) == 0
        skipped_zero = true; 
        f_0 = freq(1); 
        X_0 = X(1); 
        freq = freq(2:end); 
        X = X(2:end); 
    else
        skipped_zero = false; 
    end
    
    mask = (abs(freq) / f_rotation) .^ delta_exponent; 
    
    X_rot = mask .* X; 
    
    if skipped_zero
        freq = [f_0, freq]; 
        X_rot = [X_0, X_rot];    
    end
    
    x = real(ifft(X_rot)); 
    
    
end

% normalize
x = (x - mean(x)) / std(x); 
