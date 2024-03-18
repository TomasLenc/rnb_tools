function [amp_out] = get_amp_summary(mX, freq, frex, varargin)
% Calculate FFT magniudes sumariozed across harmonics. 
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
% method : str, optional, default='sum'
%     Method to combine the magnitudes across frex. Default is to 'sum',
%     but can also ask for 'mean'.
% 
% Returns 
% -------
% amp_sum : array_like, shape=shape(mX)(1:end-1)
%     Summed amplitude values across frex. 
% 
parser = inputParser(); 

addParameter(parser, 'method', 'sum'); % sum or mean

parse(parser, varargin{:}); 

method = parser.Results.method; 


frex_idx = ensure_row(dsearchn(ensure_col(freq), ensure_col(frex))); 

index = repmat({':'}, 1, ndims(mX)); 
index{end} = frex_idx; 

if strcmp(method, 'sum')
    amp_out = sum(mX(index{:}), ndims(mX)); 
elseif strcmp(method, 'mean')
    amp_out = mean(mX(index{:}), ndims(mX)); 
end
