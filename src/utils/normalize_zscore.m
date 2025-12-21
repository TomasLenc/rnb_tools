function z_norm = normalize_zscore(z, N, k, varargin)
% Normalize beat-related zscore to fall between -1 and 1. This is based on
% the fact that given the total number of values used in the zscoring (N),
% and the number of zscores averaged (k) to obtain the mean zscore at
% beat-related frequencies or lags, there is an analytically determined
% minimum and maximum value one can obtain. Using this min/max, we can
% normalize the zscore such that -1 reflects the smallest possible value,
% and +1 the largest possible value. This is useful to better capture the
% amount of periodization not just in relative terms (e.g. EEG vs. sound
% envlelope), but also in absolute terms (e.g. if the sound already has
% normalized zscore near 1, it means there is no room to further emphasize
% the beat periodicity - the input is already almost perfectly periodic).
% 
% Parameters
% ----------
% z : array_like 
%     Raw zscore values to normalize. 
% N : int
%     Number of values used for z-scoring. 
% k : int
%     Number of z-scores averaged to compute `z`.
% sd_method : string {'sample', 'population'}, optional, default='sample'
%     How was the SD estimated during zscoring? If it was normalized by N
%     use 'population', if by N-1 use 'sample'. 
%      
% Returns
% -------
% z_norm : float
%     Normalized zscore (range between -1 and 1). 


parser = inputParser;

addParameter(parser, 'sd_method', 'sample'); %'sample' or "population"

parse(parser, varargin{:});

sd_method = parser.Results.sd_method;

if strcmpi(sd_method, 'sample')
    % the zscores were calculated using SD normalized by N-1
    norm_factor = sqrt(  ( (N-1) * (N-k) )  / (k*N) ); 
elseif strcmpi(sd_method, 'population')
    % the zscores were calculated using SD normalized by N
    norm_factor = sqrt( (N-k) / k ); 
else 
    error('SD normalization method %s not supported', sd_method)
end

z_norm = z ./ norm_factor; 
