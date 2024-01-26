function frex = get_freq_harmonics(f0s, max_f, varargin)
% This function returns all higher harmonics (i.e. multiples) of a frequency, 
% up to the requested max frequency limit. Also, it will make sure that 
% all harmonics of any frequencies passed in the array f0s_to_exclude 
% will be excluded. 
% 
% Parameters
% ----------
% f0s : array of floats
%     Frequency(ies) to get harmonics from.
% max_f : float
%     Frequency up to which to get harmonics.
% f0s_to_exclude : array of floats
%     Frequencies for which all harmoincs overlapping with the harmonics 
%     of the requested frequncies will be excluded.
% 
% Returns
% -------
% frex : array of floats
%     All harmonics of the target frequencies, that don't overlap with the 
%     harmonics of f_harm_to_exlclude
%

parser = inputParser; 

addParameter(parser, 'f0s_to_exclude', []); 

parse(parser, varargin{:}); 

f0s_to_exclude = parser.Results.f0s_to_exclude; 


% generate all nonoveralping multiples of the requested frequencies
frex = []; 
for i_f=1:length(f0s)
    frex = [frex, f0s(i_f) : f0s(i_f) : max_f]; 
end

% make sure there is no repetition
frex = uniquetol(frex, 1e-8); 

% remove frequencies overapping with multiples of f0s_to_exclude
mask_to_rm = zeros(1, length(frex), 'like', false); 

for i_excl=1:length(f0s_to_exclude)
    
    f0_harmonics_to_exclude = [...
        f0s_to_exclude(i_excl) : f0s_to_exclude(i_excl) : max_f ...
        ]; 
        
    for li=1:length(frex)
       
        if any(abs(frex(li) - f0_harmonics_to_exclude) < 1e-8)
            mask_to_rm(li) = true; 
        end
        
    end
    
end

frex = frex(~mask_to_rm); 

