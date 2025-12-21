function [step_max_frex, slopes, frex, magns_d, magns_d_pred, magns] = ...
    get_bandwidth(mX, freq, f0, max_freq, n_harm_per_step)
% This function extracts magnitude from each harmonic at a time, and fits a
% linear model to them to estimate rate of change. 
%
% Parameters
% ----------
% header : lw6-format frequency-amplitude header
% data : lw6-format frequency-amplitude data 
%     Magnitude spectra. 
% f0 : float
%     Fundamental frequency. Harmonics of this frequency will be tested for 
%     significance. 
% max_freq : float
%     Maximum frequency taken into account. 
% n_harm_per_slope : int
%     Number of harmonics that will be successsively used to estimate the
%     slope. E.g. if f0 is equal to 1 Hz, and n_harm_per_slope is 3, the
%     function will fit a slope to the magnitude at 1, 2, 3 Hz. Then it
%     will fir a slope for 4, 5, 6 Hz. Then 7, 8, 9 Hz, etc.
%     
% Returns
% -------
% frex : cell array, shape=[n_steps, 1]
%     There is one cell per estimated slope, each containts harmonics (in
%     Hz) that were taken into account. 
% magns : cell array, shape=[n_steps, 1]
%     There is one cell per estimated slope, each containts magnitudes (per
%     harmonic) that were taken into account.
% slopes : array_like, shape=[n_steps, 1]
%     Estiamted rate of change in magnitude as a function of frequency, one
%     per each step. 

mains_frequency = 50; 

mX = ensure_row(mX); 

% get frequencies of interest
n_steps = floor(max_freq / (f0 * n_harm_per_step)); 

n_harm_all = n_harm_per_step * n_steps; 

frex_all = f0 * [1 : n_harm_all]; 

% extract magnitudes
magns_all = extract_amps(mX, freq, frex_all); 

% get derivative
magns_all_d = cumsum([magns_all]); 

% allocate maximum frequency for each channel 
frex = cell(n_steps, 1); 
magns = cell(n_steps, 1); 
magns_d = cell(n_steps, 1); 
magns_d_pred = cell(n_steps, 1); 
step_max_frex = nan(n_steps, 1); 
slopes = nan(n_steps, 1); 

for i_step=1:n_steps
    
    idx = n_harm_per_step * (i_step-1) + [1:n_harm_per_step];    ; 
    frex{i_step} = frex_all(idx); 
    magns{i_step} = magns_all(idx); 
    magns_d{i_step} = magns_all_d(idx); 
    
    step_max_frex(i_step) = frex{i_step}(end); 
    
    % remove harmomics of power line noise frequency
    if any(mod(frex{i_step}, mains_frequency) == 0)
        warning('some harmonics overalapping with 50 Hz in step %d, hope you know what you are doing...', i_step); 
    end
    
    % fit lm
    c = polyfit(frex{i_step}, magns_d{i_step}, 1);

    % extract slope
    slopes(i_step) = c(1); 
        
    % get predicted derivatives
    magns_d_pred{i_step} = c(2) + c(1) * frex{i_step}; 
    
    
end



