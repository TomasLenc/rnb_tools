function [on_off_contrast_all, on_val_all, off_val_all, on_pos_sec, off_pos_sec] = ...
                                get_on_off(x, fs, win_dur, varargin)
% Calculate the contrast in mean signal value, taken from windows starting at
% requested 'on-beat' positions and 'off-beat' positions. 
% 
% On/off-beat positions can be either defined by: 
% (1) Providing the period and phase of a regular pulse. In this case, the 
%  on-beat windows will be centered periodically starting at pulse positions. 
% The off-beat windows will be positioned at the exact midpoint between 
% on-beat positions.
% (2) Alternatively, if 'on_beat_pattern' and 'off_beat_pattern' are provided, 
% a seamlessly repeating rhythmic cycle will be assumed, and analysis windows
% will be centered on positions defined in the on- and off-beat pattern. 
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
% on_beat_pattern : array of [1, 0], optional
%     Recurring pattern of events defined on an equally-spaced periodic grid
%     of timepoints. On-beat windows will be centered on positions containing
%     '1'. 
% off_beat_pattern : array of [1, 0], optional
%     Recurring pattern of events defined on an equally-spaced periodic grid
%     of timepoints. Off-beat windows will be centered on positions containing
%     '1'. 
% grid_interval : float, optional
%     Must be provided if 'on_beat_pattern' and 'off_beat_pattern' are used.
%     Time interval in seconds between two successive points on the periodic 
%     grid (on which the patterns are defined). 
% pulse_period_sec : float, optional
%     Period of the pulse that is going to be analysed (in seconds). 
% pulse_phase_sec : float, optional, default=0
%     Phase of the pulse that is going to be analysed (in seconds). This is with
%     respect to the start of the signal (i.e. time of the first pulse position 
%     after begining of the signal).   
% on_win_offset : float, optional, default=0
%     Offset of the on-beat window in seconds, with respect to the beat
%     time. 
% off_win_offset : float, optional, default=pulse_period_sec/2
%     Offset of the off-beat window in seconds, with respect to the beat
%     time. 
% half_rectify : bool, optional, default=false
%     If true, the signal will be half-wave rectified before any metrics are
%     obtained (negative portions of the signal are set to 0). 
% full_rectify : bool, optional, default=false
%     If true, the signal will be full-wave rectified before any metrics are
%     obtained (absolute value of the signal is taken). 
% verbose : bool, optional, default=false
%     If true, diagnostics will be printed to the console. 
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

addParameter(parser, 'on_beat_pattern', []); 
addParameter(parser, 'off_beat_pattern', []); 
addParameter(parser, 'grid_interval', nan); 
addParameter(parser, 'pulse_period_sec', nan); 
addParameter(parser, 'pulse_phase_sec', 0); 
addParameter(parser, 'on_win_offset', 0); 
addParameter(parser, 'off_win_offset', nan); 
addParameter(parser, 'half_rectify', false); 
addParameter(parser, 'full_rectify', false); 
addParameter(parser, 'verbose', false); 

parse(parser, varargin{:});

on_beat_pattern = parser.Results.on_beat_pattern; 
off_beat_pattern = parser.Results.off_beat_pattern; 
grid_interval = parser.Results.grid_interval; 
pulse_period_sec = parser.Results.pulse_period_sec; 
pulse_phase_sec = parser.Results.pulse_phase_sec; 
on_win_offset = parser.Results.on_win_offset; 
off_win_offset = parser.Results.off_win_offset; 
half_rectify = parser.Results.half_rectify; 
full_rectify = parser.Results.full_rectify; 
verbose = parser.Results.verbose; 

if isnan(off_win_offset) && ~isnan(pulse_period_sec)
   off_win_offset = pulse_period_sec/2; 
end

% check we have a consistent set of parameters provided, defining either a
% method based on on/off-beat patterns (grid-based), OR a method based on
% periodic pulse positions. 
win_def_pattern_params_provided = all(...
    ~isempty(on_beat_pattern) && ...
    ~isempty(off_beat_pattern) &&...
    ~isnan(grid_interval)...
    );

win_def_period_params_provided = all(...
    ~isnan(pulse_period_sec) && ...
    ~isnan(pulse_phase_sec) &&...
    ~isnan(on_win_offset) &&...
    ~isnan(off_win_offset)...
    );

if ~xor(win_def_pattern_params_provided, win_def_period_params_provided)
    if win_def_pattern_params_provided && win_def_period_params_provided
        error(...
            'OnOff:InvalidInput', ...
            'Found some parameters for BOTH "pattern-based" and "period-based" window definition!'); 
    end
    if ~win_def_pattern_params_provided && ~win_def_period_params_provided
        error(...
            'OnOff:InvalidInput', ...
            'Did not find complete set of parameters for NEITHER "pattern-based" or "period-based" window definition!'); 
    end
end

if half_rectify && full_rectify
    error('cannot set both "half_rectify" and "full_rectify" to be true!'); 
end

win_n = round(win_dur * fs); 

% in case x is 1d column, this will fix it
x = ensure_row(x, 'verbose', false); 

shape = size(x); 

x_dur = shape(end) / fs; 

% full-wave rectify the signal if requested
if full_rectify
    x = abs(x); 
elseif half_rectify
    x = max(x, 0); 
end

% define starting times for all on- and off-beat windows 
if win_def_period_params_provided
    
    on_pos_sec = [0 + pulse_phase_sec : pulse_period_sec : x_dur] + ...
                 on_win_offset;     
    off_pos_sec = on_pos_sec + off_win_offset;
        
elseif win_def_pattern_params_provided
    
    assert(length(on_beat_pattern) == length(off_beat_pattern)); 
    
    pattern_dur = round(length(on_beat_pattern) * grid_interval, 10); 
    n_patterns_trial = floor(x_dur / pattern_dur); 
    
    on_beat_pattern_trial = repmat(on_beat_pattern, 1, n_patterns_trial);   
    on_pos_sec = grid_interval * (find(on_beat_pattern_trial) - 1); 

    off_beat_pattern_trial = repmat(off_beat_pattern, 1, n_patterns_trial);   
    off_pos_sec = grid_interval * (find(off_beat_pattern_trial) - 1); 
    
end

on_pos_sec(on_pos_sec >= x_dur) = []; 
off_pos_sec(off_pos_sec >= x_dur) = []; 
    
% go over each individual pulse position
on_val_all = nan([shape(1:end-1), length(on_pos_sec)]); 
off_val_all = nan([shape(1:end-1), length(off_pos_sec)]); 

for i_pulse=1:length(on_pos_sec)
    % get mean feature value in the on-beat window
    idx_start = round(on_pos_sec(i_pulse) * fs); 
    if idx_start + win_n > shape(end)
        break
    end
    index = repmat({':'}, 1, ndims(x)); 
    index{end} = [idx_start+1 : idx_start+win_n]; 
    val_on = mean(x(index{:}), ndims(x)); 
    % save to array
    index = repmat({':'}, 1, ndims(x)); 
    index{end} = i_pulse;     
    on_val_all(index{:}) = val_on;

end

for i_pulse=1:length(off_pos_sec)
    % get mean feature value in the off-beat window
    idx_start = round(off_pos_sec(i_pulse) * fs); 
    if idx_start + win_n > shape(end)
        break
    end
    index = repmat({':'}, 1, ndims(x)); 
    index{end} = [idx_start+1 : idx_start+win_n]; 
    val_off = mean(x(index{:}), ndims(x)); 
    % save to array
    index = repmat({':'}, 1, ndims(x)); 
    index{end} = i_pulse; 
    off_val_all(index{:}) = val_off;    
end

if verbose
    fprintf('Found %d on-beat events and %d off-beat events\n', ...
        length(~isnan(on_val_all)), length(~isnan(off_val_all))); 
end

on_val_all = nanmean(on_val_all, ndims(x));
off_val_all = nanmean(off_val_all, ndims(x)); 

% get Mitchelson contrast on-vs-off beat
on_off_contrast_all = (on_val_all - off_val_all) ./ (on_val_all + off_val_all); 










