function [s_trial, t, s_trial_dirac] = get_s(pattern, grid_ioi, fs, varargin)
% Generates whole-trial continuous signal based on a repeating pattern. 
% 
% Parameters
% ----------
% pattern : array of ones and zeros, shape=[1, N]
%     Rhythmic pattern representation with N grid-points. Events are
%     represented as 1, 0 otherwise. 
% grid_ioi : float, default=0.2
%     Time (in seconds) between successive events on the pattern grid. 
% fs : int
%     Sampling rate in Hz. 
% n_cycles : int, optional, default=1
%     Number of times the pattern should be repeated in the output signal. 
% ir : array_like, optional, default=[1]
%     Impulse response kernel. Default is Dirac impulse response. 
% jitter : float, optional, default=0
%     SD of gaussian jitter applied to the onset time of each event. 
% emph_magn : float, optional, default=0
%     Magnitude of periodic emphasis as a proportion of max amplitude. 
% emph_period : int, optional, default=4
%     Period of periodic emphasis (in units of grid points).         
% emph_phase : int, optional, default=0
%     Phase of periodic emphasis (in units of grid points).            
% 
% Returns
% -------
% s_trial : array_like
%     Time-domain representation of the resulting signal. 
% t : array_like
%     Vector of times in seconds. 
% s_trial_dirac : array_like
%     Time-domain representation using Dirac impulses, after emphasis has been
%     applied but before convolution is performed. 

parser = inputParser(); 

addParameter(parser, 'n_cycles', 1); 
addParameter(parser, 'ir', 1); 
addParameter(parser, 'emph_magn', 0); 
addParameter(parser, 'emph_period', 4); 
addParameter(parser, 'emph_phase', 0); 
addParameter(parser, 'jitter', 0); 

parse(parser, varargin{:}); 

n_cycles = parser.Results.n_cycles; 
ir = parser.Results.ir; 
emph_magn = parser.Results.emph_magn; 
emph_period = parser.Results.emph_period; 
emph_phase = parser.Results.emph_phase; 
jitter = parser.Results.jitter; 


%% 
                         
% calculate duration of one long trial in seconds
trial_dur = n_cycles * length(pattern) * grid_ioi;

% make time vector for one trial
t = [0 : 1/fs : trial_dur-1/fs];

% allocate envelope vector for the whole trial (as zeros)
s_trial_dirac = zeros(size(t));

% make ncycles copies of the rhythmic pattern
pattern_whole_trial = repmat(pattern, 1, n_cycles);

emph_magn_abs = max(abs(pattern_whole_trial)) * emph_magn; 

% go over each event in the trial
for i=1:length(pattern_whole_trial)
    
    % find the time of event onset
    event_time = (i-1)*grid_ioi;
    % apply jitter
    event_time = event_time + jitter*randn(); 
    event_time = max(event_time, 0); 
    event_time = min(event_time,trial_dur); 
    % convert to index
    event_idx = round(event_time*fs);
    
    amp = pattern_whole_trial(i); 
    % find whether there's emphasis on this grid point
    if mod((i-1-emph_phase), emph_period) == 0
        amp = amp + emph_magn_abs; 
    end

    s_trial_dirac(event_idx+1) = amp;
end

% convolve with impulse response
s_trial = conv(ir, s_trial_dirac, 'full'); 
s_trial = s_trial(1:length(t)); 