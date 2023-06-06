function x_chunked = epoch_chunks(x, fs, chunk_duration, varargin)
% Segment data into successive chunks. Assumes the first sample is at time 0.
% Time must be on the last dimension of x. 

parser = inputParser; 

addParameter(parser, 'start', 0)
addParameter(parser, 'verbose', false)

parse(parser, varargin{:})

start_time = parser.Results.start; 
verbose = parser.Results.verbose; 

shape = size(x); 

n_samples_chunk = round(chunk_duration * fs);

n_chunks = floor((shape(end) - round(start_time * fs)) / n_samples_chunk); 

x_chunked = nan([n_chunks, shape(1:end-1), n_samples_chunk]); 

t = start_time; 

for i_chunk=1:n_chunks
   
    idx_start = round(t * fs) + 1; 
    idx_end = round(t * fs) + n_samples_chunk; 

    idx_in = repmat({':'}, 1, ndims(x)); 
    idx_in{end} = [idx_start : idx_end]; 
    
    idx_out = repmat({':'}, 1, ndims(x_chunked)); 
    idx_out{1} = i_chunk; 
    
    x_chunked(idx_out{:}) = x(idx_in{:}); 
    
    t = t + chunk_duration; 
    
end

% when the input is simply a 1xN vector, it's a bit annoying to have the output
% with shape [n_chunks, 1, N]. If that's the case, let's squeeze the thing...
if ndims(x) == 2 && size(x, 1) == 1
    if verbose
        fprintf('-----------------------------------------------------------\n'); 
        fprintf('segmentation: squeezing chunked output to make life easier :)\n'); 
        fprintf('-----------------------------------------------------------\n\n'); 
    end
    x_chunked = squeeze(x_chunked); 
end

