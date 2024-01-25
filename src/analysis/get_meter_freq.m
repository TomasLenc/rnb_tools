function [freq_meter_rel, freq_meter_unrel, frex] = get_meter_freq(...
                    max_freq, varargin)

parser = inputParser; 

addParameter(parser, 'min_freq', 0); 
addParameter(parser, 'f0_to_excl', []); 

parse(parser, varargin{:}); 

min_freq = parser.Results.min_freq; 
f0_to_excl = parser.Results.f0_to_excl; 


freq_meter_rel = [1.25 : 1.25 : max_freq]; 

freq_meter_unrel = [1/2.4 : 1/2.4 : max_freq];

freq_meter_unrel = freq_meter_unrel(...
    ~ismembertol(freq_meter_unrel, freq_meter_rel, 1e-6)); 

if ~isempty(f0_to_excl)
    
    freq_meter_rel(mod(freq_meter_rel, f0_to_excl) < ...
        1e4 * eps(min(freq_meter_rel))) = []; 
    
    freq_meter_unrel(mod(freq_meter_unrel, f0_to_excl) < ...
        1e4 * eps(min(freq_meter_unrel))) = []; 

end

frex = sort([freq_meter_rel, freq_meter_unrel]);

% make sure one more time that there's no overlap between meter-rel and -unrel !!! 
assert(~any( min(abs(bsxfun(@minus, freq_meter_rel', freq_meter_unrel))) < 1e-9 ))

