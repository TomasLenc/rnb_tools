function r = correlation(x1, x2, varargin)
% Pearson correlationfor multi-dimensional arrays. Time must be last dimension!

parser = inputParser(); 

addParameter(parser, 'method', 'pearson'); 

parse(parser, varargin{:}); 

method = parser.Results.method; 

%%

if strcmpi(method, 'pearson')
    
    nd = ndims(x1); 

    x1_norm = x1 - mean(x1, nd); 
    x2_norm = x2 - mean(x2, nd); 

    num = sum(x1_norm .* x2_norm, nd); 

    den = sqrt(sum(x1_norm.^2, nd) .* sum(x2_norm.^2, nd)); 

    r = num ./ den; 
    
elseif strcmpi(method, 'spearman')
        
    shape = size(x1); 
    r = nan(shape(1:end-1), 1); 

    % prepare index vector as cell, e.g. {1,1,1,1,':'} 
    nv = ndims(x1) - 1;  % exclude last dimension
    idx_while_loop = [repmat({1}, 1, nv), {':'}]; 
    max_size_per_dim = size(x1); % size of each dimension
    ready = false; 

    while ~ready
        % apply function 
        r(idx_while_loop{:}) = corr(...
            squeeze(x1(idx_while_loop{:}))',...
            squeeze(x2(idx_while_loop{:}))', ...
            'type', 'Spearman' ...
            ); 
        % Update the index vector:
        % Assume that the WHILE loop is ready
        ready = true;       
         % Loop through dimensions
        for k = 1:nv        
            % Increase current index by 1
            idx_while_loop{k} = idx_while_loop{k} + 1;   
            % Does it exceed the length of current dim?
            if idx_while_loop{k} <= max_size_per_dim(k) 
               % No, WHILE loop is not ready now
               ready = false;  
               % idx_while_loop(k) increased successfully, leave "for k" loop
               break;         
            end
            % Reset idx_while_loop{k}, proceed to next k
            idx_while_loop{k} = 1;          
        end
    end
       
else
    error('method "%s" not implemented', method); 
end