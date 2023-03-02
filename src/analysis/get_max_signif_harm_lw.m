function [max_signif_freq, n_chan_signif_last_iter] = get_max_signif_harm_lw(...
                                header, data, f0, n_harm_per_test, snr_bins, ...
                                varargin)
% This function keeps testing for significance of higher and higher groups of 
% harmonics. It returns the max harmonic (in Hz) above which the mean harmonic 
% is not significantly above the noise floor anymore. 
%
% Parameters
% ----------
% header : lw6-format frequency-amplitude header
% data : lw6-format frequency-amplitude data 
%     Magnitude spectra. 
% f0 : float
%     Fundamental frequency. Harmonics of this frequency will be tested for 
%     significance. 
% n_harm_per_test : int
%     Number of harmonics that will be tested for each iteration. E.g. if f0 is 
%     equal to 1 Hz, and n_harm_per_test is 3, the function will start testing 
%     whether the average response at 1, 2, 3 Hz is above noise. Then it will 
%     test 4, 5, 6 Hz. Then 7, 8, 9 Hz, etc.
% snr_bins : [int, int]
%     Min and max frequency bin from each side around the target frequency that 
%     will be used for noise estimation.
% thr_p : float, optional, default=0.01
%     Pvalue threshold for testing significance of signal > noise. 
% n_successive_nonsignif_steps_to_stop : int, optional, default=1
%     Number of successive iterations that must be nonsignificant, in order to 
%     stop testing the particular electrode. 
% mains_frequency : float, optional, default=60
%     Power line noise frequency. (in Mexico it's 60 Hz). 
%     
% Returns
% -------
% max_signif_harm : array of floats, shape=[channel, 1]
%     Maximum frequency up to which the response is significant, separately for
%     each channel.

parser = inputParser; 

addParameter(parser, 'thr_p', 0.01); 
addParameter(parser, 'mains_frequency', 60); 
addParameter(parser, 'n_successive_nonsignif_steps_to_stop', 1); 

parse(parser, varargin{:});

thr_p = parser.Results.thr_p; 
mains_frequency = parser.Results.mains_frequency; 
n_successive_nonsignif_steps_to_stop = parser.Results.n_successive_nonsignif_steps_to_stop; 

% calculate threshold zscore
thr_z = norminv(1 - thr_p);

% prepare array of frequencies
freq = [0 : header.datasize(end) - 1] * header.xstep; 

% allocate counter of groups of harmonics we've tested
c = 0;
% allocate maximum frequency for each channel 
max_signif_freq = zeros(header.datasize(2), 1);
max_reached = false(header.datasize(2), 1);
n_successive_nonsignificant_steps = zeros(header.datasize(2), 1);

% check that we don't have multiple trials
if header.datasize(1) > 1
    error('We can only process 1 average trial, you passed %d', header.datasize(1));
end

while true
    
    frex = f0 * n_harm_per_test * c + [1 : n_harm_per_test] * f0;
    
    % stop if we're at the end of our frequency vector and still significant...
    if frex(end) > freq(end - snr_bins(2))
        break
    end
    
    % remove harmomics of power line noise frequency
    frex(mod(frex, mains_frequency) == 0) = [];
    
    % get zscore snr for each channel
    z_snr = get_z_snr(data, freq, frex, ...
                      snr_bins(1),...
                      snr_bins(2)); 
    
    z_snr = ensure_col(z_snr);
    
    signif_chan_mask = z_snr > thr_z;
    
    % reset counter for channels that are significant in this iteration, and 
    % are still in the game
    n_successive_nonsignificant_steps(~max_reached & signif_chan_mask) = 0;
    
    % increate the counter of channels that are still in the game and were
    % nonsignificant in this iteration
    n_successive_nonsignificant_steps(~max_reached) = ...
                        n_successive_nonsignificant_steps(~max_reached) + ...
                        ~signif_chan_mask(~max_reached); 
    
    % take out of the game channels that haven't been significnat in 
    % n_nonsignif_steps_before_stop succesive steps, so that we don't update them anymore 
    max_reached(n_successive_nonsignificant_steps >= n_successive_nonsignif_steps_to_stop) = true;
    
    % if all channels are nonsignificant by now, exit the looop
    if all(max_reached)
        break 
    end
    
    c = c+1;
    
    % for significant channels that we're still updating, write the current
    % maximum frequency
    max_signif_freq(~max_reached & signif_chan_mask) = f0 * n_harm_per_test * c;
    
end

assert(all(n_successive_nonsignificant_steps(max_reached) == ...
            n_successive_nonsignif_steps_to_stop))

assert(all(n_successive_nonsignificant_steps(~max_reached) < ...
            n_successive_nonsignif_steps_to_stop))
