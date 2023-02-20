function [header, data] = make_dataset(varargin)
% This function makes an artificial continuous dataset for testing. 
% 
% Parameters
% ----------
% fs : int
%     Sampling rate.
% trial_dur : float
%     Duration of a single trial. 
% event_sequence : cell of strings (Optional)
%     If passed, the dataset will contain the exact sequence of
%     event codes marking successie trials. If absent, you have to specify
%     n_cond and n_trials instead. 
% n_cond : int (Optional)
%     Number of simulated conditions (each will have an event code starting 
%     from '1', '2', ...) 
% n_trials : int (Optional)
%     Number of trials per condition. Trials will be shuffled.
% pause_dur : float or array([float, float])
%     Number of seconds between successive trials. If scalar, the pause
%     duration will be fixed. If array of two scalars, the pause duration will
%     be randomly sampled from the given range. 
% n_chans : int 
%     Number of channels.
% start_recording_buffer : float 
%     Pause before the first trial in seconds.
% end_recording_buffer : float 
%     Pause before the last trial in seconds.
% 
% Returns
% -------
% header : letswave6 header
% data : letswave6 data
%

parser = inputParser;

addParameter(parser, 'fs', 512);
addParameter(parser, 'trial_dur', 33.6);
addParameter(parser, 'event_sequence', {});
addParameter(parser, 'n_cond', 10);
addParameter(parser, 'n_trials', 10);
addParameter(parser, 'pause_dur', [2, 5]);
addParameter(parser, 'n_chans', 64);
addParameter(parser, 'start_recording_buffer', 10);
addParameter(parser, 'end_recording_buffer', 10);

parse(parser, varargin{:});

fs = parser.Results.fs;
trial_dur = parser.Results.trial_dur;
n_trials = parser.Results.n_trials;
event_sequence = parser.Results.event_sequence;
n_cond = parser.Results.n_cond;
n_chans = parser.Results.n_chans;
pause_dur_params = parser.Results.pause_dur;
start_recording_buffer = parser.Results.start_recording_buffer;
end_recording_buffer = parser.Results.end_recording_buffer;

%% prepare data

trial_codes = {};
if ~isempty(event_sequence)
    trial_codes = event_sequence;
else
    for i_cond=1:n_cond
        trial_codes = [trial_codes, repmat({sprintf('%d', i_cond)}, 1, n_trials)];
    end
    trial_codes = randsample(trial_codes, length(trial_codes)); 
end

if numel(pause_dur_params) == 1
    pause_durs = repmat(pause_dur_params, 1, length(trial_codes));    
elseif numel(pause_dur_params) == 2
    pause_durs = rand(1, length(trial_codes) - 1) * ...
                     abs(diff(pause_dur_params)) + ...
                     min(pause_dur_params);
else
    error('pause duration must be either 1 or 2 element vector!');
end

total_dur = start_recording_buffer + ...
            end_recording_buffer +...
            sum(pause_durs(1:end-1)) + ...
            length(trial_codes) * trial_dur;

N = round(total_dur * fs); 
data = 0.1 * rand([1, n_chans, 1, 1, 1, N]); 

t = start_recording_buffer;
c = 1;

for i_trial=1:length(trial_codes)
        
    if i_trial == length(trial_codes)
        pause_dur = 0;
    else
        pause_dur = pause_durs(i_trial);
    end
    
    idx = round(t * fs) + 1;
    data(:, :, :, :, :, idx) = str2num(trial_codes{i_trial});
    
    events(c).code = trial_codes{i_trial};
    events(c).latency = t;
    events(c).epoch = 1;
    
    t = t + pause_dur + trial_dur;
    c = c + 1;

end


%% get header

header = make_header('data', data, ...
                     'x_start', 0, ...
                     'fs', fs, ...
                     'n_chans', n_chans, ...
                     'events', events);









