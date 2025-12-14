function [flux, fs_flux] = get_spectral_flux(s, fs, varargin)

parser = inputParser; 

addParameter(parser, 'fs_target', []); 
addParameter(parser, 'N_target', []); 

parse(parser, varargin{:}); 

fs_target = parser.Results.fs_target; 
N_target = parser.Results.N_target; 

% ----------------------------------------------------------------------
% MIR toolbox 
% ----------------------------------------------------------------------

% addpath(genpath('~/projects_git/mirtoolbox')); 
% 
% 
% fname = '/Users/tomaslenc/projects_git/Intracerebral/Intracerebral_NaturalMusic_stimuli/Bugz.wav'; 
% 
% audio = miraudio(fname); 
% 
% % Compute spectral flux between consecutive frames
% % 50 ms window, 25 ms hop
% flux = mirflux(audio);
% 
% flux_vals = mirgetdata(flux);
% dt_flux = 0.025; 
% fs_flux = 1 / dt_flux; 
% t_flux = dt_flux + [0 : length(flux_vals)-1] * dt_flux; 
% 
% % resample 
% % flux = interp1(t_flux, flux_vals, t, 'linear');
% [p, q] = rat(fs / fs_flux);   
% pad = zeros(1, round(dt_flux * fs_flux)); 
% flux = resample([pad, flux_vals, pad], p, q);
% 
% figure
% plot(t_flux, flux_vals)
% hold on 
% plot(t, flux)

%%

% ----------------------------------------------------------------------
% Based on  Section 6.1.2 (spectral-based novelty detection) of book
% Fundamentals of Music Processing by Meinard Müller.
% ----------------------------------------------------------------------

% Steps:
% 1. Read audio
% 2. STFT
% 3. Magnitude ? Log compression
% 4. Positive temporal differences ? Spectral Flux
% 5. Local-average subtraction ? Enhanced Spectral Flux (optional)


if size(s, 2) == 2
    s = mean(s, 2);   % Convert to mono if stereo
end

win_length = 2048; % round(0.050 * fs_audio); 
hop_size = 512; % round(win_length / 2); 
gamma = 1;           % Log compression factor (? ? 1)
M = 10;          % Local averaging half-window (in frames)

% Zero-pad the beginning to allow centered windowing
s = [zeros(win_length,1); s; zeros(win_length,1)];

% Number of frames
n_frames = floor((length(s) - win_length)/hop_size) + 1;

% Hann window
w = hann(win_length, 'periodic');

% Preallocate: K = N/2
K = floor(win_length/2);
X = zeros(K+1, n_frames);

ptr = 1;
for n = 1:n_frames
    frame = s(ptr : ptr+win_length-1) .* w;
    X_frame = fft(frame);
    X(:,n) = abs(X_frame(1:K+1));   % keep only 0...Nyquist
    ptr = ptr + hop_size;
end

% compression
Y = log(1 + gamma * X);

% Temporal derivative
dY = diff(Y,1,2);   % Y(:,n+1) - Y(:,n)

% Keep only positive increases
dY(dY < 0) = 0;

% Flux(n) = sum over frequency bins
flux = sum(dY, 1);

% Flux time index aligns with frame 2..numFrames
t_flux = ((0:n_frames-1) * hop_size) / fs;
t_flux = t_flux(2:end);
t_flux = [0, t_flux]; 

flux = [0, flux]; 

% resample and cut 
fs_flux = 1 / (hop_size / fs); 
if ~isempty(fs_target)
    [p, q] = rat(fs_target / fs_flux);   
    flux = resample(flux, p, q);
    fs_flux = fs_target; 
end
if ~isempty(N_target)
    flux = flux(1:N_target); 
end

% mu = movmean(flux, 2*M+1);     % local average
% flux_enh = max(flux - mu, 0);
% 
