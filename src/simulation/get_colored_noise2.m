function x = get_colored_noise2(shape, fs, exponent)
% Generates multidimensional colored noise. 
% 
% Parameters
% ----------
% shape : array_like
%     Dimensions of the output (time is always the last dimension!). 
% fs : int
%     Sampling rate.
% exponent : float
%     Exponent f^beta (careful! not 1/(f^beta)) so don't mess up the sign ;)
% 
% Returns
% -------
% x : array_like
%     Colored noise with the requested shape and exponent. 

% in order to call bsxfun we need time on the first axis...we'll flip it all
% back at the end! 
shape = flip(shape); 

N = shape(1); 
T = N/fs; 
freq = [0 : N-1] / T; 

%% method from Cohen 

hN = ceil(N / 2) - 1; 

% random magnitudes decreasing with 1/f shape 
magnitudes = bsxfun(@times, rand([hN, shape(2:end)]), freq(2:hN+1)' .^ (exponent/2));

% uniform random phase for each frequency 
phases = exp(1j * 2 * pi * rand([hN, shape(2:end)])); 

% prepare zeros that will go to 0 and hN frequency bin
zeros_to_fill = zeros([1, shape(2:end)]); 

% get flipped magnidues and phases
magnitudes_rev = flip(magnitudes, 1); 
phases_rev = conj(flip(phases, 1)); 

% create whole complex spectrum
if mod(N, 2)==0    
    
    full_magnitudes = cat(1, ...
                          zeros_to_fill, magnitudes, ...
                          zeros_to_fill, magnitudes_rev); 
                      
    full_phases = cat(1,...
                      zeros_to_fill, phases, ...
                      zeros_to_fill, phases_rev);     
                  
else
    
    full_magnitudes = cat(1, ...
                          zeros_to_fill, magnitudes, ...
                          magnitudes_rev); 
                      
    full_phases = cat(1,...
                      zeros_to_fill, phases, ...
                      phases_rev);     
                      
end

X = full_magnitudes .* full_phases; 

x = real(ifft(X, [], 1)); 

x = (x - mean(x, 1)) ./ std(x, [], 1); 

x = permute(x, [length(shape):-1:1]); 


%% method from neurodsp (not implemented for multi-D yet)

% parser = inputParser; 
% addParameter(parser, 'f_rotation', 1); 
% parse(parser, varargin{:});
% f_rotation = parser.Results.f_rotation; 
% 
% x = randn(1, N); 
% 
% 
% % Delta exponent is divided by two, as the FFT output is in units of amplitude not power
% delta_exponent = exponent / 2; 
% 
% X = fft(x); 
% 
% if freq(1) == 0
%     skipped_zero = true; 
%     f_0 = freq(1); 
%     X_0 = X(1); 
%     freq = freq(2:end); 
%     X = X(2:end); 
% else
%     skipped_zero = false; 
% end
% 
% mask = (abs(freq) / f_rotation) .^ delta_exponent; 
% 
% X_rot = mask .* X; 
% 
% if skipped_zero
%     freq = [f_0, freq]; 
%     X_rot = [X_0, X_rot];    
% end
% 
% x = real(ifft(X_rot)); 
% 
% x = (x - mean(x)) / std(x); 
% 

