function x_str = vec_to_str(vec, varargin)

parser = inputParser(); 

addParameter(parser, 'format', '%.0f'); 
addParameter(parser, 'sep', ''); 

parse(parser, varargin{:}); 

format_str = parser.Results.format; 
sep_str = parser.Results.sep; 

vec = ensure_row(vec); 

x_str = strjoin(cellfun(@(x) num2str(x, format_str), ...
            num2cell(vec, 1), 'uni', 0), sep_str); 