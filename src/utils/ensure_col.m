function x = ensure_col(x, varargin)

parser = inputParser; 

addParameter(parser, 'verbose', false); 

parse(parser, varargin{:});

verbose = parser.Results.verbose; 


if isempty(x)
    return
elseif iscolumn(x)
    return
elseif isrow(x)
    x = x.';  
    return
else
    if verbose
        warning('too many dimensions... keeping as it is...'); 
    end
end
