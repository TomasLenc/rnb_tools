function patterns = generate_rhythms(N, K, varargin)
% GENERATE_RHYTHMS
%
% Generates rhythmic patterns with:
%   - length = N
%   - exactly K events (1s)
%   - no groups of events larger than max_group (default = 3)
%   - no groups of silences larger than max_silence (default = Inf)
%   - first element = 1
%   - optional removal of cyclic rotations
%
% USAGE:
%   patterns = generate_rhythms(N, K)
%   patterns = generate_rhythms(N, K, 'max_group', 2)
%   patterns = generate_rhythms(N, K, 'max_silence', 2)
%   patterns = generate_rhythms(N, K, 'remove_rotations', false)
%
% INPUTS:
%   N : integer
%       pattern length
%
%   K : integer
%       number of events
%
% OPTIONAL PARAMETERS:
%   'max_group'        : integer (default = 3)
%       maximum allowed size of consecutive-event groups
%
%   'max_silence'      : integer (default = Inf)
%       maximum allowed size of consecutive-silence groups
%
%   'remove_rotations' : logical (default = true)
%       if true, removes cyclically equivalent patterns
%
% OUTPUT:
%   patterns : matrix (num_patterns x N)

%% Input parsing

p = inputParser;

addRequired(p, 'N', @(x) isnumeric(x) && isscalar(x) && x > 0);
addRequired(p, 'K', @(x) isnumeric(x) && isscalar(x) && x >= 0);

addParameter(p, 'max_group', 3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'max_silence', Inf, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'remove_rotations', true, @(x) islogical(x) && isscalar(x));

parse(p, N, K, varargin{:});

max_group        = p.Results.max_group;
max_silence      = p.Results.max_silence;
do_remove_rot    = p.Results.remove_rotations;

%% Generate combinations

comb = nchoosek(1:N, K);
numComb = size(comb, 1);

valid_patterns = [];

fprintf('going over %d patterns... \n', numComb); 
for i = 1:numComb
    ptn = zeros(1, N);
    ptn(comb(i,:)) = 1;

    if ptn(1) ~= 1
        continue;
    end

    if has_large_group(ptn, max_group, 1)
        continue;
    end

    if has_large_group(ptn, max_silence, 0)
        continue;
    end

    valid_patterns = [valid_patterns; ptn];
end

%% Remove cyclic duplicates (optional)

if do_remove_rot
    patterns = remove_rotations(valid_patterns);
else
    patterns = valid_patterns;
end

fprintf('%d patterns passed the criteria... \n', size(patterns, 1)); 

end


%% Helper functions

function flag = has_large_group(p, max_len, value)
% Returns true if pattern has a group of 'value' larger than max_len (circular)

if isinf(max_len)
    flag = false;
    return;
end

pp = [p p];
group_len = 0;
flag = false;

for i = 1:length(pp)
    if pp(i) == value
        group_len = group_len + 1;
        if group_len > max_len
            flag = true;
            return;
        end
    else
        group_len = 0;
    end
end

end


function unique_patterns = remove_rotations(patterns)
% Remove cyclically equivalent patterns

numP = size(patterns,1);
keep = true(numP,1);

for i = 1:numP
    if ~keep(i)
        continue;
    end

    p = patterns(i,:);

    rotations = zeros(length(p), length(p));
    for r = 0:length(p)-1
        rotations(r+1,:) = circshift(p, [0 r]);
    end

    for j = i+1:numP
        if keep(j) && any(ismember(rotations, patterns(j,:), 'rows'))
            keep(j) = false;
        end
    end
end

unique_patterns = patterns(keep,:);

end