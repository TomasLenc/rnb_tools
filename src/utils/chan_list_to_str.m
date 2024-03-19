function out_str = chan_list_to_str(labels, varargin)

parser = inputParser(); 

addParameter(parser, 'sep', ', '); 

parse(parser, varargin{:}); 

sep = parser.Results.sep; 

if length(labels) == 1
    out_str = labels{1}; 
else
    out_str = strjoin(labels, sep); 
end