function amps = extract_amps(mX, freq, frex, varargin)
% Extract values from FFT. 
% 
% Parameters
% ----------
% mX : array_like, shape=[..., frequency]
%     Noise-subtracted magnitude spectra with frequency as the last
%     dimension.
% freq : array_like
%     Frequencies for the FFT. 
% frex : array_like
%     Frequencies to take the snippets around. 
% method : str, optional, default='nearest'
%     Method to extract the value. 'nearest' looks at the value in a bin
%     closest to the requested frequency. 'max' extracts maximum value
%     within the closest bin and `n_bins` on each side. 'mean' takes the
%     average of all those bins. The latter two options can be useful if
%     leakage is expected.
% n_bins : int
%     Number of bins on each side of the bin closest to the requested
%     frequency will be considered. 
% 
% Returns 
% -------
% amp_sum : array_like, shape=shape(mX)(1:end-1)
%     Summed amplitude values across frex. 
% 
parser = inputParser(); 

addParameter(parser, 'method', 'nearest'); % nearest, max, mean
addParameter(parser, 'n_bins', 0); 

parse(parser, varargin{:}); 

method = parser.Results.method; 
n_bins = parser.Results.n_bins; 

%%

frex_idx = ensure_row(dsearchn(ensure_col(freq), ensure_col(frex))); 

if strcmp(method, 'nearest')
    
    index = repmat({':'}, 1, ndims(mX)); 
    index{end} = frex_idx; 
    amps = mX(index{:}); 

else

    tmp_size = size(mX); 
    tmp_size = [tmp_size(1:end-1), length(frex), 2*n_bins + 1]; 
    tmp = nan(tmp_size); 

    for i_f=1:length(frex)

        index_mX = repmat({':'}, 1, ndims(mX)); 
        index_mX{end} = [frex_idx(i_f)-n_bins : frex_idx(i_f)+n_bins]; 

        index_tmp = repmat({':'}, 1, ndims(tmp)); 
        index_tmp{end-1} = i_f; 

        tmp(index_tmp{:}) = mX(index_mX{:}); 

    end

    if strcmp(method, 'max')
        amps = max(tmp, [], ndims(tmp)); 
    elseif strcmp(method, 'mean')
        amps = mean(tmp, ndims(tmp)); 
    else 
        error('method not implemented', method); 
    end

end

    
    
    
    
    
    
    
    
    
    
    
    