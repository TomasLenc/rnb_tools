function [on_off_contrast_all, on_val_all, off_val_all] = get_on_off(...
                                x, fs, pulse_period_sec, win_dur, varargin...
                                )
% Calculate phase stability of narrowband-filtered signal phase across pulse 
% positions. 
% 
% Parameters
% ----------
% x : array_like, shape=[..., time]
%     Time-domain signal to analyse. Time must be the last dimension. 
% fs : float
%     Sampling rate (samples/s)
% pulse_period_sec : float
%     Period of the pulse that is going to be analysed (in seconds). 
% win_dur : float
%     Duration of the analysis window starting at 'on'- and 'off'-beat
%     positions. Mean value of the signal within these windows will be taken to
%     calculate the on/off contrast. 
% on_win_offset : float, optional, default=0
%     Offset of the on-beat window in seconds, with respect to the beat
%     time. 
% off_win_offset : float, optional, default=pulse_period_sec/2
%     Offset of the off-beat window in seconds, with respect to the beat
%     time. 
% pulse_phase_sec : float, optional, default=0
%     Phase of the pulse that is going to be analysed (in seconds). This is with
%     respect to the start of the signal (i.e. time of the first pulse position 
%     after begining of the signal).   
% 
% Returns 
% -------
% on_off_contrast_all : float
%     Michelson contrast of mean signal value in the on- vs. off-beat windows. 
% on_val_all : float
%     Mean value of the signal averaged across all on-beat windows. 
% off_val_all : float
%     Mean value of the signal averaged across all off-beat windows. 

parser = inputParser; 

addParameter(parser, 'pulse_phase_sec', 0); 
addParameter(parser, 'on_win_offset', 0); 
addParameter(parser, 'off_win_offset', pulse_period_sec/2); 

parse(parser, varargin{:});

pulse_phase_sec = parser.Results.pulse_phase_sec; 
on_win_offset = parser.Results.on_win_offset; 
off_win_offset = parser.Results.off_win_offset; 

%% 

win_n = round(win_dur * fs); 

shape = size(x); 

x_dur = size(x, ndims(x)) / fs; 
pulse_pos_sec = [0 + pulse_phase_sec : pulse_period_sec : x_dur-1/fs]; 

on_val_all = nan([shape(1:end-1), length(pulse_pos_sec)]); 
off_val_all = nan([shape(1:end-1), length(pulse_pos_sec)]); 
on_off_contrast_all = nan([shape(1:end-1), length(pulse_pos_sec)]); 

% go over each individual pulse position
for i_pulse=1:length(pulse_pos_sec)

    % get mean feature value in the on-beat window
    idx_start = round((pulse_pos_sec(i_pulse) + on_win_offset) * fs); 

    index = repmat({':'}, 1, ndims(x)); 
    index{end} = [idx_start+1 : idx_start+win_n]; 
    val_on = mean(x(index{:}), ndims(x)); 

    % get mean feature value in the off-beat window
    idx_start = round((pulse_pos_sec(i_pulse) + off_win_offset) * fs); 

    index = repmat({':'}, 1, ndims(x)); 
    index{end} = [idx_start+1 : idx_start+win_n]; 
    val_off = mean(x(index{:}), ndims(x)); 

    % get Mitchelson contrast on-vs-off beat
    index = repmat({':'}, 1, ndims(x)); 
    index{end} = i_pulse; 
    
    on_val_all(index{:}) = val_on;
    off_val_all(index{:}) = val_off;
    on_off_contrast_all(index{:}) = (val_on - val_off) ./ (val_on + val_off); 
    
end

on_val_all = mean(on_val_all, ndims(x)); 
off_val_all = mean(off_val_all, ndims(x)); 
on_off_contrast_all = mean(on_off_contrast_all, ndims(x)); 










