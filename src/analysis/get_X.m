function [X, freq, fres] = get_X(x, fs, varargin)
% Compute Complex DFT of a time-domain signal. Time must be the last
% dimension. 
parser = inputParser(); 

addParameter(parser, 'dc_zero', false); % set DC (0 Hz) to zero

parse(parser, varargin{:}); 

dc_zero = parser.Results.dc_zero; 


if iscell(x)
    
    assert(length(unique(cellfun(@(x) length(x), x, 'uni', 1))) == 1)
    N = size(x{1}, ndims(x{1})); 
    hN = floor(N/2); 
    
    X = cellfun(@(x) fft(x, [], ndims(x)), x, 'uni', 0); 
    
    index = repmat({':'}, 1, ndims(x{1})); 
    index{end} = [1 : hN]; 
    X = cellfun(@(X) X(index{:}), X, 'uni', 0); 
    
    if dc_zero
        for i=1:length(X)
            index = repmat({':'}, 1, ndims(X{i})); 
            index{end} = 1; 
            X{i}(index{:}) = 0; 
        end
    end
        
else
    
    N = size(x, ndims(x)); 
    hN = floor(N/2); 
    X = fft(x, [], ndims(x)); 
    
    index = repmat({':'}, 1, ndims(x)); 
    index{end} = [1 : hN]; 
    X = X(index{:}); 

    if dc_zero
        index = repmat({':'}, 1, ndims(X)); 
        index{end} = 1; 
        X(index{:}) = 0; 
    end

end

freq = [0 : hN-1] / N * fs; 

fres = 1 / N * fs; 