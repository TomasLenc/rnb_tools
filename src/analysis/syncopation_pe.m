function [C_min, best_phase, accents, C_all] = syncopation_pe(pattern, period, varargin)
% SYNCOPATION_PE Compute Povel & Essens (1985) syncopation (C-score)
%
%   [C_min, best_phase, accents, C_all] = syncopation_pe(pattern, period)
%   computes the counterevidence score (C-score) for a binary rhythmic
%   pattern according to Povel & Essens (1985).
%
%   The algorithm evaluates all possible phase alignments of a metrical
%   grid
%   and returns the minimal counterevidence.
%
%   C = weight * O + U
%
%   where:
%       O = number of omitted accents (strong beats with silence)
%       U = number of unaccented onsets (events without phenomenal accent)
%
%   INPUTS
%   -------
%   pattern : vector of 0/1
%       Binary rhythmic pattern (1 = onset, 0 = silence)
%
%   period : integer
%       Metrical period (group size), e.g., 4
%
%   OPTIONAL PARAMETERS
%   --------------------
%   'weight'   : scalar (default = 4)
%       Weight for omitted accents relative to unaccented onsets
%
%   'circular' : logical (default = true)
%       If true, the pattern is treated as circular.
%       If false, zero-padding is used for edge handling.
%
%   OUTPUTS
%   --------
%   C_min     : scalar
%       Minimum counterevidence score across all phase alignments
%
%   best_phase : integer
%       Phase (0-based) yielding minimum C-score
%
%   accents   : vector
%       Binary vector marking phenomenal accents
%
%   C_all : vector
%       C-score for each phase (length = period)
%       Index corresponds to phase (0-based via index-1)
%
%   NOTES
%   -----
%   Phenomenal accents are defined as:
%       1) A note followed by silence
%       2) The start of a group of >= 3 consecutive notes
%
%   The function evaluates all phase shifts of the metric grid and returns
%   the best-fitting alignment (minimum C).
%
%   REFERENCE
%   ---------
%   Povel, D.-J., & Essens, P. (1985). "Perception of temporal patterns."
%   Music Perception, 2(4), 411?440.

%% --- Parse inputs ---

p = inputParser;
addParameter(p, 'weight', 4);
addParameter(p, 'circular', true);
parse(p, varargin{:});

W = p.Results.weight;
circular = p.Results.circular;

pattern = pattern(:)' > 0;
N = length(pattern);

%% --- Padding ---
% Handles edge events

if circular
    padded = [pattern pattern pattern];
else
    padded = [zeros(1,N) pattern zeros(1,N)];
end


%% --- Accent detection ---

accents = zeros(1, N);

for i = 1:N
    idx = i + N;

    if padded(idx) == 1

        % Rule 1: followed by silence
        if padded(idx + 1) == 0
            accents(i) = 1;
            continue;
        end

        % Rule 2: start of >=3-note group
        if padded(idx - 1) == 0 && ...
           padded(idx + 1) == 1 && ...
           padded(idx + 2) == 1
            accents(i) = 1;
        end
    end
end

%% --- C-score for each phase ---

C_all = zeros(1, period);

for phase = 0:(period - 1)

    idx = (phase + 1):period:N;

    O = sum(pattern(idx) == 0);
    U = sum(pattern(idx) == 1 & accents(idx) == 0);

    C_all(phase + 1) = W * O + U;
end

%% --- Best phase ---

[C_min, best_idx] = min(C_all);
best_phase = best_idx - 1;


